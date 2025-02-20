import subprocess

if __name__ == "__main__":
    cmd = "poetry run rich-click --output svg afragmenter --help > help_message.svg"

    subprocess.run(cmd, shell=True)

    with open("help_message.svg", 'r') as svg:
        content = svg.read()
        content = content.replace(">Rich</text>", ">AFragmenter</text>")

    with open("help_message.svg", 'w') as svg:
        svg.write(content)
    
    print("Help message saved as help_message.svg")