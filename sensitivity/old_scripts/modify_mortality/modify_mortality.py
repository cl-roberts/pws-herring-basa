import csv as csv
import pandas as pd
def modify_mortality(new_value):
    # Extract tabular data from the CTL file
    data = []
    with open('../model/PWS_ASA(par).ctl', 'r') as file:
        for line in file:
            if not line.startswith('#') and line.strip():  # Skip comments or empty lines
                columns = line.strip().split()  # Assuming space-separated values
                data.append(columns)

    # Write to a CSV file#
        with open('output.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(data)


    # Load CSV file into DataFrame
    df = pd.read_csv('output.csv', nrows = 31, skiprows = [0])
    df.loc[-1] = df.columns  # Add the column names as the first row
    df.index = df.index + 1  # Shift the index to move the row to the top
    df = df.sort_index().reset_index(drop=True)  # Sort and reset the index
    df.columns = ['init_value', 'lower_bound', 'upper_bound', 'est_phz', 'prior_type', "p1", "p2", 
            'fun_type', '#', 'number', 'parameter']

    index_keep = df.index[df['parameter'] == 'Z_0_8'].tolist()
    df.iloc[index_keep, 0] = new_value

    return df
        