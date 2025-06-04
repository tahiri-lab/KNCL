#!/usr/bin/env python

"""
Script Name: <overlap_tree_pipeline.py>
Last Modified: 2025-04-20
Description: This Python script automates the process of creating datasets of biologically meaningful partially overlapping phylogenetic trees with branch lengths.
"""

import pandas as pd
import random
import os
import sys
import math
import time
import shutil
import zipfile
import requests
from Bio import Phylo

# For Selenium
from selenium import webdriver
from selenium.common.exceptions import TimeoutException, NoSuchElementException, WebDriverException
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options

# Parse command-line arguments
if len(sys.argv) != 5:
    print("Usage: python overlap_tree_pipeline.py <species_group> <n> <number_of_trees> <email>")
    sys.exit(1)

species_group = sys.argv[1].lower()
n = int(sys.argv[2])
number_of_trees = int(sys.argv[3])
email = sys.argv[4]

allowed_species_groups = ['amphibians', 'birds', 'mammals', 'sharks', 'squamates']

if species_group not in allowed_species_groups:
    print(f"Species group must be one of: {', '.join(allowed_species_groups)}")
    sys.exit(1)

if number_of_trees < 10 or number_of_trees > 1000 or number_of_trees % 10 != 0:
    print("Error: The number of trees must be between 10 and 1000 and divisible by 10.")
    sys.exit(1)

t = number_of_trees // 10  # Number of trees to select per subset

if not os.path.exists('all_species_lists.csv'):
    print("Error: 'all_species_lists.csv' file not found. Please ensure the file is present in the current directory.")
    sys.exit(1)

def calculate_n_for_k(k, p_values=[0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]):
    n = k  # Start with k species in the first subset
    previous_new_species = 0  # Track the new species in the previous step
    for p in p_values:
        # Calculate the number of common species
        common_species = (2 * k * p) / (1 + p)
        rounded_common_species = round(common_species)
        new_species = k - rounded_common_species
        actual_new_species = new_species - previous_new_species
        n += math.ceil(actual_new_species)
        previous_new_species = new_species
    return int(n)

def find_k_from_n(n, max_k=1000, p_values=[0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]):
    for k in range(1, max_k + 1):
        calculated_n = calculate_n_for_k(k, p_values)
        if calculated_n >= n:
            return k
    return None

# Alternative:
# p_min = min(p_values)
# k = round(n * (1 + p_min) / 2)

# Load the data from the CSV file
data = pd.read_csv('all_species_lists.csv')

# Extract the taxa for each group
amphibians = data['Amphibians'].dropna().tolist()
birds = data['Birds'].dropna().tolist()
mammals = data['Mammals'].dropna().tolist()
sharks = data['Sharks'].dropna().tolist()
squamates = data['Squamates'].dropna().tolist()

species_dict = {
    'amphibians': amphibians,
    'birds': birds,
    'mammals': mammals,
    'sharks': sharks,
    'squamates': squamates
}

species_list = species_dict[species_group]

if len(species_list) < n:
    print(f"Error: Not enough species in {species_group} to select {n} species.")
    sys.exit(1)

# Randomly select n species
final_species = random.sample(species_list, n)

# Create a DataFrame
result_df = pd.DataFrame({species_group.capitalize(): final_species})

# Save to 'selected_species.csv'
output_file = 'selected_species.csv'
result_df.to_csv(output_file, index=False)

def create_overlapping_subsets(df, p_values=[0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]):
    for group in df.columns:
        species_list = df[group].dropna().tolist()
        n = len(species_list)

        # Calculate the subset size k
        k = find_k_from_n(n)
        if k is None:
            print(f"Could not find a valid k for group {group}")
            continue

        # Create subsets
        subsets = {f"Subset {i}": [] for i in range(1, 11)}
        start_index = 0  # Start index remains fixed

        # Start with the first subset
        subsets['Subset 1'] = species_list[start_index:start_index + k]
        current_position = start_index + k  # Update the current position to move down the list

        for i, p in enumerate(p_values, start=2):
            # Calculate common and new species
            common_species = round((2 * k * p) / (1 + p))
            new_species = k - common_species

            # Select k species starting from start_index + new_species
            subset_start = start_index + int(new_species)
            subset_end = subset_start + k

            # Handle wrapping around the species list
            subset_species = species_list[subset_start:subset_end]
            if len(subset_species) < k:
                subset_species += species_list[:k - len(subset_species)]
            # alternative: subset_species = [species_list[(j % n)] for j in range(subset_start, subset_start + k)]

            subsets[f"Subset {i}"] = subset_species

            # Update the current position for the next subset
            current_position = subset_start + k

        # Pad all subsets to ensure they are the same length
        max_len = len(species_list)
        for key in subsets:
            subsets[key] += [None] * (max_len - len(subsets[key]))  # Pad with None values if necessary

        # Save as CSV
        output_df = pd.DataFrame(subsets)
        output_file = f"{group}_overlapping_subsets.csv"
        output_df.to_csv(output_file, index=False)
        print(f"Saved file for {group}: {output_file}")

# Load your CSV file
df = pd.read_csv('selected_species.csv')

# Run the script to create overlapping subsets
create_overlapping_subsets(df)

# Mapping group names to the value attributes used on the webpage
GROUP_MAPPING = {
    "amphibians": "amphibiantree",
    "birds": "birdtree",
    "mammals": "mammaltree",
    "sharks": "sharktree",
    "squamates": "squamatetree"
}

def setup_driver():
    chrome_options = Options()
    chrome_options.add_argument("--headless=new")  # Use newer headless mode
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--disable-gpu")
    chrome_options.add_argument("--window-size=1920x1080")
    driver = webdriver.Chrome(options=chrome_options)
    return driver

def clean_folder(folder_name):
    """
    Delete all contents in the folder if it exists.
    """
    if os.path.exists(folder_name):
        print(f"Cleaning folder: {folder_name}")
        shutil.rmtree(folder_name)  # Remove the folder and its contents
        os.makedirs(folder_name)  # Recreate an empty folder
    else:
        os.makedirs(folder_name)  # Create the folder if it does not exist

def submit_tree_request(driver, species_list, email, group):
    """
    Submit a tree request and capture job ID from the webpage.
    """
    driver.get('https://vertlife.org/phylosubsets/')

    try:
        # Handle radio button selection for species group
        if group.lower() != "amphibians":  # Amphibians is selected by default
            print(f"Selecting species group: {group}")
            group_value = GROUP_MAPPING.get(group.lower())
            if group_value:
                group_radio_button = WebDriverWait(driver, 30).until(
                    EC.element_to_be_clickable((By.XPATH, f"//input[@type='radio' and @value='{group_value}']"))
                )
                group_radio_button.click()
            else:
                print(f"Group {group} not found in the mapping.")
        else:
            print("Skipping species group selection since Amphibians is the default.")

        # Handle any potential pop-ups or modals here

        # Paste species into the textarea
        species_textarea = WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.ID, 'selected'))
        )
        species_textarea.clear()
        species_textarea.send_keys('\n'.join(species_list))  # Use newline separator for species list

        print(f"Email being used: {email}")
        # Enter email address
        email_input = driver.find_element(By.ID, 'email')
        email_input.clear()
        email_input.send_keys(email)

        print("Submitting tree request...")
        # Click 'Get Trees'
        get_trees_button = driver.find_element(By.ID, 'btnGetTrees')
        get_trees_button.click()

        # Wait for job ID
        print("Waiting for Job ID...")
        job_id_element = WebDriverWait(driver, 60).until(
            EC.presence_of_element_located((By.XPATH, "//div[@id='status']//strong"))
        )
        job_id = job_id_element.text.strip()  # Extract the job ID
        print(f"Job ID: {job_id}")
        return job_id

    except Exception as e:
        print(f"Error during tree request: {str(e)}")
        return None

def download_zipfile_using_job_id(job_id, folder_name, max_retries=20):
    """
    Directly construct the URL from the job ID and download the file, retry up to max_retries if it fails.
    """
    download_url = f"https://data.vertlife.org/pruned_treesets/{job_id}/{job_id}.zip"
    print(f"Downloading from: {download_url}")

    retries = 0
    while retries < max_retries:
        response = requests.get(download_url)

        if response.status_code == 200:
            zip_filename = os.path.join(folder_name, f"{job_id}.zip")

            # Save the zip file
            with open(zip_filename, 'wb') as f:
                f.write(response.content)

            try:
                with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
                    print(f"Extracting files from {zip_filename}...")
                    for file in zip_ref.namelist():
                        # Rename extracted files using job_id to ensure uniqueness
                        file_extension = os.path.splitext(file)[1]  # Get the file extension (e.g., .nex)
                        new_file_name = f"{job_id}{file_extension}"  # Construct new file name using job_id
                        target_path = os.path.join(folder_name, new_file_name)

                        # Extract and rename the file
                        zip_ref.extract(file, folder_name)
                        os.rename(os.path.join(folder_name, file), target_path)
                        print(f"Extracted: {target_path}")
            except zipfile.BadZipFile:
                print(f"Error: {zip_filename} is not a valid zip file.")
            return True
        else:
            retries += 1
            print(f"Failed to download from {download_url}. Status code: {response.status_code}. Retrying... ({retries}/{max_retries})")
            time.sleep(30)
    print(f"Failed to download file after {max_retries} attempts.")
    return False

def automate_process(input_file, email):
    # Extract group name dynamically from the input file name
    group_name = os.path.basename(input_file).split('_')[0]
    print(f"Group Name: {group_name}")

    # Create and clean folder for the group
    folder_name = f"{group_name}_nexus"
    clean_folder(folder_name)  # Clean the folder before starting the process

    # Load the input data (species subsets)
    df = pd.read_csv(input_file)

    job_ids = {}

    # Step 1: Request trees and collect job IDs
    for subset_name in df.columns:
        success = False
        retries = 0
        while not success and retries < 15:
            try:
                print(f"Requesting trees for {subset_name}")
                species_list = df[subset_name].dropna().tolist()
                # Initialize the driver for each request
                driver = setup_driver()
                job_id = submit_tree_request(driver, species_list, email, group_name)
                driver.quit()  # Quit the driver after the request
                if job_id:
                    job_ids[subset_name] = job_id
                    success = True
                else:
                    print(f"Retrying {subset_name} due to missing job ID.")
                    retries += 1
                    time.sleep(20)
            except Exception as e:
                print(f"Exception during processing subset {subset_name}: {e}")
                retries += 1
                time.sleep(20)

        # Wait a bit before the next request
        time.sleep(10)

    # Added Delay Before Downloading Files
    print("Waiting for 60 seconds before starting downloads to ensure files are ready...")
    time.sleep(60)  # Wait for 60 seconds

    # Step 2: Download the results using the job ID
    for subset_name, job_id in job_ids.items():
        print(f"Processing download for {subset_name}")
        download_success = download_zipfile_using_job_id(job_id, folder_name)
        if not download_success:
            print(f"Failed to download trees for {subset_name}")

input_file = f"{species_group}_overlapping_subsets.csv"
automate_process(input_file, email)

def convert_nexus_to_newick(nexus_file, t):
    """
    Convert trees from a Nexus file to Newick format, randomly select t trees.
    """
    # Parse the Nexus file
    try:
        trees = list(Phylo.parse(nexus_file, "nexus"))
    except Exception as e:
        print(f"Error parsing {nexus_file}: {e}")
        return []

    # Check if we have enough trees in the file
    if len(trees) < t:
        print(f"Warning: {nexus_file} has fewer than {t} trees. Selecting all available trees.")
        t = len(trees)

    # Randomly select t trees
    selected_trees = random.sample(trees, t)

    # Convert the selected trees to Newick format and clean up ':0.00000' at the end
    selected_newick_trees = []
    for tree in selected_trees:
        newick_str = tree.format('newick').strip()
        # Remove ':0.00000' or any ':<number>' before the final ';'
        if newick_str.endswith(":0.00000;"):
            newick_str = newick_str[:-9] + ';'
        selected_newick_trees.append(newick_str)

    return selected_newick_trees

def process_nexus_files(input_dir, output_file, t):
    """
    Process all Nexus files in the input directory, convert trees to Newick,
    randomly select t trees from each file, and combine them into one output file.
    """
    all_selected_trees = []

    # Loop through all .nex files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".nex"):
            nexus_file = os.path.join(input_dir, filename)
            selected_trees = convert_nexus_to_newick(nexus_file, t)
            all_selected_trees.extend(selected_trees)

    # Save all selected trees into one file without blank lines
    with open(output_file, 'w') as out_file:
        for tree in all_selected_trees:
            if tree:  # Ensure no blank lines
                out_file.write(tree + '\n')

    print(f"Saved {len(all_selected_trees)} Newick trees to {output_file}")

input_dir = f"./{species_group}_nexus/"  # Adjust the path as needed
output_file = f"overlapping_dataset_{species_group}.txt"  # Path to the final combined output file

process_nexus_files(input_dir, output_file, t)

# Citations
citations = {
    'amphibians': """Amphibians: Jetz, W., & Pyron, R. A. (2018). The interplay of past diversification and evolutionary isolation with present imperilment across the amphibian tree of life. Nature ecology & evolution, 2(5), 850-858.""",
    'birds': """Birds: Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). The global diversity of birds in space and time. Nature, 491(7424), 444-448.""",
    'mammals': """Mammals: Upham, N. S., Esselstyn, J. A., & Jetz, W. (2019). Inferring the mammal tree: species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLoS biology, 17(12), e3000494.""",
    'sharks': """Sharks: Stein, R. W., Mull, C. G., Kuhn, T. S., Aschliman, N. C., Davidson, L. N., Joy, J. B., ... & Mooers, A. O. (2018). Global priorities for conserving the evolutionary history of sharks, rays and chimaeras. Nature ecology & evolution, 2(2), 288-298.""",
    'squamates': """Squamates: Tonini, J. F. R., Beard, K. H., Ferreira, R. B., Jetz, W., & Pyron, R. A. (2016). Fully-sampled phylogenies of squamates reveal evolutionary patterns in threat status. Biological Conservation, 204, 23-31."""
}

print(f"Dataset saved to {output_file}")
print("Please cite the following source for the data used:")
print(citations[species_group])
print("Website with comprehensive data: https://vertlife.org/data/")
