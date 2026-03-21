import pandas as pd

# Step 1: Input the number of files to read
n = int(input("Enter the number of files to read: "))

# Step 2: Initialize an empty DataFrame to store data from each file
combined_data = pd.DataFrame()

# Step 3: Loop through each file, read the data, and add it as a new column
for i in range(1, n + 1):
    txt_file_path = f'ChromosomeSpacerFreq{i}.txt'  # Construct the filename
    #txt_file_path = f'ChromosomeSpacerFreqPromoter{i}.txt'  # Construct the filename
    try:
        # Read each .txt file
        data = pd.read_csv(txt_file_path, delimiter='\t', header=None)
        
        # Add this data as a new column to the combined DataFrame
        combined_data[f'ChromosomeSpacerFreq{i}'] = data[0]  # Assumes each file has a single column
        #combined_data[f'ChromosomeSpacerFreqPromoter{i}'] = data[0]  # Assumes each file has a single column

        print(f"Successfully added {txt_file_path} as a column.")
    except FileNotFoundError:
        print(f"Warning: {txt_file_path} not found. Skipping this file.")
    except Exception as e:
        print(f"Error reading {txt_file_path}: {e}")

# Step 4: Save the combined data to an .xlsx file
xlsx_file_path = 'combined_output.xlsx'
combined_data.to_excel(xlsx_file_path, index=False)

for i in range(1, n + 1):
   # txt_file_path = f'ChromosomeSpacerFreq{i}.txt'  # Construct the filename
    txt_file_path = f'ChromosomeSpacerFreqPromoter{i}.txt'  # Construct the filename
    try:
        # Read each .txt file
        data = pd.read_csv(txt_file_path, delimiter='\t', header=None)
        
        # Add this data as a new column to the combined DataFrame
       # combined_data[f'ChromosomeSpacerFreq{i}'] = data[0]  # Assumes each file has a single column
        combined_data[f'ChromosomeSpacerFreqPromoter{i}'] = data[0]  # Assumes each file has a single column

        print(f"Successfully added {txt_file_path} as a column.")
    except FileNotFoundError:
        print(f"Warning: {txt_file_path} not found. Skipping this file.")
    except Exception as e:
        print(f"Error reading {txt_file_path}: {e}")

# Step 4: Save the combined data to an .xlsx file
xlsx_file_path = 'combined_output2.xlsx'
combined_data.to_excel(xlsx_file_path, index=False)

print(f"Data from {n} files successfully written to {xlsx_file_path}")
