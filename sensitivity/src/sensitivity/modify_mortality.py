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
    f = open(dir_model + '/PWS_ASA(par).ctl', 'r')
    lines = f.readlines()
    # write new line
    lines[17] = "    " + '%.2f' % new_value + "    0.05    1.50     -2      0      0     0     0   # 3:  Z_0_8\n"
    f.close()

    # write ctl
    f = open(dir_model + '/PWS_ASA(par).ctl', 'w')
    f.writelines(lines)
    f.close()
