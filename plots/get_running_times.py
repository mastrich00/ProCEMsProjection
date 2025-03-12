from datetime import datetime

# Define the two datetime strings
datetime_str1 = "2025-02-18 11:17:42.307299"
datetime_str2 = "2025-02-18 11:26:04.723575"

# Convert the strings to datetime objects
datetime_obj1 = datetime.strptime(datetime_str1, "%Y-%m-%d %H:%M:%S.%f")
datetime_obj2 = datetime.strptime(datetime_str2, "%Y-%m-%d %H:%M:%S.%f")

# Calculate the difference
time_difference = datetime_obj2 - datetime_obj1
print(time_difference.total_seconds())
