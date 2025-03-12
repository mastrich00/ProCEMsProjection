import csv
import glob

# Define the target actions in the second column that trigger a replacement
# target_actions = {"Redund proCEMs", "Enumerate projected Cone", "Postprocessing proCEMs"}
target_actions = {"Enumerate projected Cone"}
# Find all CSV files in the current directory
for csv_file in glob.glob("testResults/measurements/*_ecoli_*/*csv.csv"):
    # Read the file contents
    with open(csv_file, 'r', newline='') as f:
        reader = csv.reader(f)
        rows = list(reader)
    
    # Process each row
    for row in rows:
        if len(row) >= 2 and row[1] in target_actions:
            row[1] = "Enumerating ProCEMs"
    
    # Write the modified data back to the same file
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(rows)