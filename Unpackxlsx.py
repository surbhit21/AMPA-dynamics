import os
import pandas as pd

root_folder = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Published/K_willig/"

excel_file = os.path.join(root_folder,'TableS4_relto_Figure4.xlsx')
all_sheets = pd.read_excel(excel_file, sheet_name=None)
sheets = all_sheets.keys()
snames = []
for sheet_name in sheets:
    sheet = pd.read_excel(excel_file, sheet_name=sheet_name)
    snames.append("Fig4_{}".format(sheet_name))
    sheet.to_csv(os.path.join(root_folder,"Fig4/Fig4_{}.csv".format(sheet_name)), index=False)

print(snames)