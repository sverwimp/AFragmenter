import numpy as np
from scipy.optimize import linear_sum_assignment
import re


def get_max_residue_index(domain_strings):
    """
    Determines the maximum residue index across all provided domain strings.
    Handles negative starting positions by treating them as 0.
    """
    max_idx = -1
    for domain_string in domain_strings:
        if not domain_string:
            continue
        individual_domains_str = domain_string.split(',')
        for dom_str in individual_domains_str:
            ranges_str = dom_str.split('_')
            for r_str in ranges_str:
                if '-' in r_str:
                    try:
                        start, end = parse_range_with_negatives(r_str)
                        max_idx = max(max_idx, end)
                    except ValueError:
                        continue
    return max_idx


def parse_domain_string_np(domain_string, max_residue_index):
    """
    Parses a domain string into a list of 1D NumPy boolean arrays (masks).
    Negative starting positions are treated as 0.
    """
    domains = []
    if not domain_string or domain_string.strip() == "":
        return domains
    
    # Split by '|' to get individual domains
    individual_domains_str = domain_string.split(',')
    
    for dom_str in individual_domains_str:
        if not dom_str.strip():  # Skip empty domain strings
            continue
            
        # Initialize a boolean array for the current domain
        current_domain_mask = np.full(max_residue_index + 1, False, dtype=bool)
        
        # Split by '_' to get fragmented segments of the same domain
        ranges_str = dom_str.split('_')
        for r_str in ranges_str:
            if not r_str.strip():  # Skip empty range strings
                continue
            if '-' in r_str:  # Ensure it's a range
                try:
                    # Handle negative numbers correctly
                    start, end = parse_range_with_negatives(r_str)
                    if start <= end and end >= 0 and start <= max_residue_index:
                        # Clamp indices to valid range
                        start_idx = max(0, start)
                        end_idx = min(max_residue_index, end)
                        if start_idx <= end_idx:
                            current_domain_mask[start_idx : end_idx + 1] = True
                except ValueError:
                    print(f"Warning: Could not parse range '{r_str}' in domain string '{domain_string}'")
                    continue
        
        # Only add if the domain mask actually contains True values
        if np.any(current_domain_mask):
            domains.append(current_domain_mask)
    return domains


def parse_range_with_negatives(range_str):
    """
    Parse a range string that may contain negative numbers using regex.
    Negative start values are set to 0.
    """
    # Pattern to match: optional minus, digits, dash, optional minus, digits
    pattern = r'^(-?\d+)-(-?\d+)$'
    match = re.match(pattern, range_str)
    
    if not match:
        raise ValueError(f"Invalid range format: {range_str}")
    
    start = int(match.group(1))
    end = int(match.group(2))
    
    # Set negative start values to 0
    if start < 0:
        start = 0

    if end < 0:
        end = 0
    
    return start, end


def calculate_iou_np(mask1, mask2):
    """Calculates IoU for two NumPy boolean masks."""
    # If both masks are entirely False (empty domains), consider IoU as 1.0 (perfect match of nothing)
    # This might be context-dependent. For domain prediction, usually domains are non-empty.
    if not np.any(mask1) and not np.any(mask2):
        return 1.0
    
    # Intersection: True where both masks are True
    intersection_mask = np.logical_and(mask1, mask2)
    intersection_size = np.sum(intersection_mask)
    
    # Union: True where either mask is True
    union_mask = np.logical_or(mask1, mask2)
    union_size = np.sum(union_mask)
    
    if union_size == 0:
        return 0.0 # No common residues, no total span
    
    return intersection_size / union_size


def iou_chain(true_domain_string, predicted_domain_string):
    """
    Calculates the IoU chain score for protein domain prediction based on the given formula,
    using optimal bipartite matching (Hungarian algorithm) for assigning predicted domains
    to true domains, respecting unique assignments.
    Uses NumPy boolean arrays for domain representation and IoU calculation.

    Args:
        true_domain_string (str): A string representing the true domains.
                                  Example: "0-76_185-227,238-244"
        predicted_domain_string (str): A string representing the predicted domains.
                                       Example: "0-81_183-244,82-180"

    Returns:
        tuple: (IoU chain score, list of IoU values per true domain)
    """
    # Determine the global maximum residue index for consistent array sizing
    max_idx = get_max_residue_index([true_domain_string, predicted_domain_string])
    # If max_idx is -1, it means no valid ranges were found in either string.
    # In such a case, the sequence length for masks would be 0, leading to issues.
    # We set a minimum size of 1 if no residues are found, though this indicates empty inputs.
    if max_idx < 0:
        max_idx = 0 # Ensures masks are at least size 1 for empty cases
        
    print(f"true_domain_string: {true_domain_string}")
    print(f"predicted_domain_string: {predicted_domain_string}")

    T_masks = parse_domain_string_np(true_domain_string, max_idx)
    P_masks = parse_domain_string_np(predicted_domain_string, max_idx)

    n_dom_T = len(T_masks)
    n_dom_P = len(P_masks)

    if n_dom_T == 0:
        return 0.0, [] # No true domains, score is 0.

    # Calculate total true domain size. This is the denominator for the weighting term.
    total_true_domain_size = sum(np.sum(t_mask) for t_mask in T_masks)
    if total_true_domain_size == 0:
        # This occurs if all 'true' domains are empty (e.g., ",," or empty ranges).
        # In typical domain prediction, true domains are expected to have size.
        return 0.0, [0.0] * n_dom_T

    # If there are no predicted domains, no matching can occur, so sum of terms will be 0.
    if n_dom_P == 0:
        return 0.0, [0.0] * n_dom_T

    # Create a cost matrix for bipartite matching.
    # Dimensions: (number of true domains) x (number of predicted domains)
    # We want to maximize IoU, so we define cost as (1.0 - IoU).
    cost_matrix = np.zeros((n_dom_T, n_dom_P), dtype=float) 

    for i in range(n_dom_T):
        for j in range(n_dom_P):
            # Calculate IoU for each possible true-predicted domain pair
            iou_val = calculate_iou_np(T_masks[i], P_masks[j])
            cost_matrix[i, j] = 1.0 - iou_val # Store as cost

    # Apply the Hungarian algorithm (linear_sum_assignment)
    # It returns two arrays: row_ind (indices of rows selected) and col_ind (indices of columns selected).
    # These represent the optimal 1-to-1 pairings that minimize the total cost.
    matched_true_indices, matched_pred_indices = linear_sum_assignment(cost_matrix)

    # Initialize a list to store the IoU value for each true domain.
    # By default, unmatched true domains will have an IoU of 0 for their term.
    iou_per_true_domain = [0.0] * n_dom_T

    # Populate iou_per_true_domain with the IoU values from the optimal matching.
    for i_true_dom, i_pred_dom in zip(matched_true_indices, matched_pred_indices):
        # Retrieve the IoU value for the matched pair (T_masks[i_true_dom], P_masks[i_pred_dom])
        # We can either recalculate it or get it from the cost_matrix (1 - cost).
        iou_val = 1.0 - cost_matrix[i_true_dom, i_pred_dom]
        iou_per_true_domain[i_true_dom] = iou_val

    # Sum the terms according to the IoU Chain formula.
    sum_terms = 0.0
    for i in range(n_dom_T):
        t_mask = T_masks[i]
        
        # The IoU term for Ti comes from the calculated matched_iou_per_true_domain.
        iou_term = iou_per_true_domain[i]
        
        # The weighting term: |Ti| / sum(|Tj|)
        # np.sum(t_mask) gives the number of 'True' residues (the size of the domain).
        weight_term = np.sum(t_mask) / total_true_domain_size
        
        sum_terms += (iou_term * weight_term)

    return sum_terms, iou_per_true_domain


if __name__ == "__main__":
    # Example usage
    true_domain = "0-76_185-227,238-244"
    predicted_domain = "0-81_183-244,82-180"
    
    iou_chain_score, iou_per_domain = iou_chain(true_domain, predicted_domain)
    print(f"IoU Chain Score: {iou_chain_score:.4f}")

    high_iou_domains = [iou for iou in iou_per_domain if iou >= 0.8]
    print(f"Number of true domains with IoU >= 0.8: {len(high_iou_domains)}")
    print(f"IoU values per true domain: {[f'{iou:.4f}' for iou in iou_per_domain]}")