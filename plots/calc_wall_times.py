import glob
import csv
from datetime import datetime

def calculate_duration(csv_file):
    """Calculate time difference between Start and End records in a CSV file"""
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        rows = list(reader)
        
        if len(rows) < 2:
            return None  # Not enough data
        
        # Extract start and end times
        start_row = rows[0]
        end_row = rows[-1]
        
        if start_row[0] != 'Start' or end_row[0] != 'End':
            return None  # Invalid file format
            
        try:
            start_time = datetime.fromisoformat(start_row[3])
            end_time = datetime.fromisoformat(end_row[3])
            return (end_time - start_time).total_seconds()
        except (IndexError, ValueError) as e:
            print(f"Error processing {csv_file}: {str(e)}")
            return None

# Find all CSV files recursively in measurements directory
csv_files = glob.glob('testResults/measurements/**/*.csv', recursive=True)

# Process files and calculate durations
durations = {}
for file_path in csv_files:
    duration = calculate_duration(file_path)
    if duration is not None:
        durations[file_path] = duration

# Print results
print("File Durations Analysis:")
for file_path, duration in durations.items():
    print(f"- {file_path}:")
    print(f"  Total duration: {duration:.2f} seconds")
    print(f"               ({duration/60:.2f} minutes)")
    print(f"               ({duration/3600:.2f} hours)")

# Optional: Add to existing plotting code
# You can use these durations in your plot titles or annotations