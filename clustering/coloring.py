import time
from pathlib import Path

def prepare_VMD_colors(
        init_lines: list[str],
        colors_input: Path = "./colors.idx",
        output: Path = "./output.pdb") -> None:
    with open(colors_input) as incolors:
        color_lines = [ic.strip() for ic in incolors.readlines()]

    with open(output, "w") as tofile:

        for idx, line in enumerate(init_lines):
            # print(f"[INIT|] {line}")

            if idx != 0:
                try:
                    color_type = int(color_lines[idx-1]) + 2 # -1 because we need to skip the first line and +2 because it started from 0
                except IndexError:
                    # print("Default color.")
                    color_type = 1

                new_line = line[:24] + f"{color_type:>2}" + line[27:]
            else:
                new_line = line
            # print(f"[FINAL] {new_line}")
            tofile.write(f"{new_line}\n")

            # print("-------------")
            # time.sleep(1)

    print(f"[{output}] is created.")

if __name__ == "__main__":
    with open("./input.pdb") as infile: init_lines = [il.strip() for il in infile.readlines()]
    prepare_VMD_colors(init_lines)

    print("Done.")
