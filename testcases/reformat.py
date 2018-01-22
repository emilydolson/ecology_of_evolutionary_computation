import sys

result = []

with open(sys.argv[1], "r") as filename:
    file_lines = filename.readlines()

    for line in file_lines:
        full_line = line.split(",")
        full_line[0] = full_line[0].strip("[]")
        inputs = full_line[0].split()
        inputs.append(full_line[1].strip())
        result.append(",".join(inputs))


with open(sys.argv[1], "w") as filename:
    filename.write("\n".join(result))