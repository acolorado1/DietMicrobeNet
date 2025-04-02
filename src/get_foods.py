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

# Display the DataFrame with the selected organisms
if not kegg_df.empty:
    st.write("### Selected Organism DataFrame", kegg_df)
else:
    st.write("No organisms selected.")

# Convert the DataFrame to CSV for download
csv = kegg_df.to_csv(index=False)

# Provide the download button
if not kegg_df.empty:
    st.download_button(
        label="Download Selected Organisms DataFrame as CSV",
        data=csv,
        file_name="kegg_organisms_dataframe.csv",
        mime="text/csv"
    )

# start the FooDB selection
st.header("FooDB Food Item Selection")

######################################### THIS IS A LOCAL DIRECTORY
food_df = pd.read_csv('/Users/burkhang/Code_Projs/DietMicrobeNet/Data/food.csv')

# subset to columns needed 
food_df = food_df[['id', 'name', 'name_scientific']]

# Allow the user to select multiple foods
selected_foods = st.multiselect("Select Food Item", food_df['name'].tolist())

# Create an empty DataFrame to hold the selected foods' data
foodb_df = pd.DataFrame(columns=["id", "name", "name_scientific"])

# Filter the DataFrame to include only selected foods
foodb_df = food_df[food_df["name"].isin(selected_foods)]

# Display the DataFrame with the selected organisms
if not foodb_df.empty:
    st.write("### Selected Food Item DataFrame", foodb_df)
else:
    st.write("No foods selected.")

# Convert the DataFrame to CSV for download
csv = foodb_df.to_csv(index=False)

# Provide the download button
if not foodb_df.empty:
    st.download_button(
        label="Download Selected Food Items DataFrame as CSV",
        data=csv,
        file_name="foodb_foods_dataframe.csv",
        mime="text/csv"
    )

# put a kill button so user can close application easily 
st.write("""
**Once you're done hit the *Stop Server* button wait 5 seconds close the tab**
""")
if st.button("Stop Server"):
    st.warning("Stopping the Streamlit server...")
    os.kill(os.getpid(), signal.SIGTERM)