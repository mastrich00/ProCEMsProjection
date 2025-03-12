import csv
import glob

# Define the action to remove
target_action = "Redund proCEMs"

# Find all CSV files in the directory
for csv_file in glob.glob("testResults/measurements/*_excel_*/*times.csv"):
    # Read and filter the file contents
    with open(csv_file, 'r', newline='') as f:
        reader = csv.reader(f)
        filtered_rows = [row for row in reader if not (len(row) >= 2 and row[1] == target_action)]
    
    # Write the filtered data back to the same file
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(filtered_rows)