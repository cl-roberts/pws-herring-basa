"""
Modify assumed natural mortality in BASA
"""

def modify_mortality(new_value, dir_model):
    """
    This function takes a positive number as input and modifies a BASA input file
    so that it changes the assumed rate of natural mortality.

    :param float new_value: New value for natural mortality
    :param str dir_model: Path to BASA model directory
    """

    if new_value < 0:
        raise ValueError("Instantaneous natural mortality cannot be negative. Change new_value")

    # read .ctl
    with open(dir_model + '/PWS_ASA(par).ctl', 'r', encoding="utf-8") as data:
        lines = data.readlines()
        # write new line
        new_line = "    " + f'{new_value:.2f}' + "    0.05    1.50     -2      0"
        new_line = new_line + "      0     0     0   # 3:  Z_0_8\n"
        lines[17] = new_line
        data.close()

    # write ctl
    with open(dir_model + '/PWS_ASA(par).ctl', 'w', encoding="utf-8") as data:
        data.writelines(lines)
        data.close()
