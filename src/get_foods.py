import streamlit as st
import requests
import os
import signal
import pandas as pd

# Fetch KEGG organism data
def fetch_kegg_organisms():
    """
    Fetch the list of KEGG organisms.
    Returns:
        List of tuples containing (code, scientific_name, common_name).
    """
    url = "http://rest.kegg.jp/list/organism"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.text.splitlines()
        organisms = []
        for line in data:
            parts = line.split("\t")
            if len(parts) == 4:
                t_id, code, scientific_name, taxonomy = parts
                organisms.append((code, scientific_name, taxonomy))
        return organisms
    else:
        st.error("Failed to fetch KEGG organism data.")
        return []

# Streamlit app
st.title("Whole Genome and/or Metabolome Representation of Diet")

st.write("""
This app allows you to customize a set of food items using KEGG organisms or FooDB (food metabolomes) that can represent a patients diet. 
Dataframes will be created based on user selections which can be downloaded as a CSV file. 
         
Disclaimer: While KEGG is kept up-to-date, version 1.0 of FooDB is used from the CSV file added April7, 2020.
""")

st.header("KEGG Organism Selection")

st.write("To ease search, if using common name, surround it with parentheses, e.g. *(cow)*")
# Fetch the list of KEGG organisms
organisms = fetch_kegg_organisms()

# Create a DataFrame from the organism list for easy lookup
organism_df = pd.DataFrame(organisms, columns=["Code", "ScientificName", "Taxonomy"])

# Allow the user to select multiple organisms
selected_organisms = st.multiselect("Select KEGG Organisms", organism_df["ScientificName"].tolist())

# Create an empty DataFrame to hold the selected organisms' data
kegg_df = pd.DataFrame(columns=["Code", "ScientificName", "Taxonomy"])

# Filter the DataFrame based on selected organism names
kegg_df = organism_df[organism_df["ScientificName"].isin(selected_organisms)]

# Display the DataFrame with the selected organisms and allow user input
if not kegg_df.empty:
    st.write("### Enter values (between 1 and 100) for each selected organism")

    new_column_values = []
    for index, row in kegg_df.iterrows():
        value = st.number_input(
            label=f"Enter value for {row['ScientificName']} ({row['Code']})",
            min_value=1,
            max_value=100,
            step=1,
            key=f"input_{index}"
        )
        new_column_values.append(value)

    # Add the new user input column
    kegg_df["food_frequency"] = new_column_values

    # Display updated DataFrame
    st.write("### Selected Organism DataFrame with Input", kegg_df)

    # Enable download
    csv = kegg_df.to_csv(index=False)
    st.download_button(
        label="Download Selected Organisms DataFrame as CSV",
        data=csv,
        file_name="kegg_organisms_dataframe.csv",
        mime="text/csv", 
        key="download_keeg"
    )
else:
    st.write("No organisms selected.")

# start the FooDB selection
st.header("FooDB Food Item Selection")

# Load the food CSV (adjust the path if needed) ####################Local Directory 
#food_df = pd.read_csv('/Users/burkhang/Code_Projs/DietMicrobeNet/Data/food.csv')

# Get the current directory of the script
script_dir = os.path.dirname(__file__)
file_path = os.path.join(script_dir, '..', 'data', 'food.csv')

# Load CSV
food_df = pd.read_csv(file_path)

# Subset to relevant columns
food_df = food_df[['id', 'name', 'name_scientific']]

# Allow the user to select multiple foods
selected_foods = st.multiselect("Select Food Item", food_df['name'].tolist())

# Filter the DataFrame to include only selected foods
foodb_df = food_df[food_df["name"].isin(selected_foods)]

# Display the DataFrame with user input column
if not foodb_df.empty:
    st.write("### Enter values (between 1 and 100) for each selected food")

    food_values = []
    for index, row in foodb_df.iterrows():
        value = st.number_input(
            label=f"Enter value for {row['name']} ({row['id']})",
            min_value=1,
            max_value=100,
            step=1,
            key=f"food_input_{index}"
        )
        food_values.append(value)

    foodb_df["food_frequency"] = food_values

    st.write("### Selected Food Item DataFrame with Input", foodb_df)

    csv = foodb_df.to_csv(index=False)

    st.download_button(
        label="Download Selected Food Items DataFrame as CSV",
        data=csv,
        file_name="foodb_foods_dataframe.csv",
        mime="text/csv",
        key="download_foodb"
    )
else:
    st.write("No foods selected.")

# put a kill button so user can close application easily 
st.write("""
**Once you're done hit the *Stop Server* button wait 5 seconds close the tab**
""")
if st.button("Stop Server"):
    st.warning("Stopping the Streamlit server...")
    os.kill(os.getpid(), signal.SIGTERM)