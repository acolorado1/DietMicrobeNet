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
st.title("KEGG Organism Selection and Downloadable DataFrame")

st.write("""
This app allows you to select a KEGG organism and add it to a DataFrame. 
You can then download the updated DataFrame as a CSV file.
""")

# Fetch the list of KEGG organisms
organisms = fetch_kegg_organisms()

# Create a list of organism names for user selection
organism_names = [f"{org[1]} ({org[0]})" for org in organisms]  # Display name and code

# Allow the user to select multiple organisms
selected_organisms = st.multiselect("Select KEGG Organisms", organism_names)

# Create an empty DataFrame to hold the selected organisms' data
df = pd.DataFrame(columns=["Code", "ScientificName", "Taxonomy"])

# Add the selected organisms to the DataFrame
for selected_organism in selected_organisms:
    selected_organism_code = selected_organism.split('(')[-1].strip(')')
    
    # Find the organism details from the full list of organisms
    selected_details = next((org for org in organisms if org[0] == selected_organism_code), None)
    
    if selected_details:
        # Add the organism details to the DataFrame
        df = df._append({
            "Code": selected_details[0],
            "ScientificName": selected_details[1],
            "Taxonomy": selected_details[2]
        }, ignore_index=True)

# Display the DataFrame with the selected organisms
if not df.empty:
    st.write("### Selected Organism DataFrame", df)
else:
    st.write("No organisms selected.")

# Convert the DataFrame to CSV for download
csv = df.to_csv(index=False)

# Provide the download button
if not df.empty:
    st.download_button(
        label="Download Selected Organisms DataFrame as CSV",
        data=csv,
        file_name="kegg_organisms_dataframe.csv",
        mime="text/csv"
    )


# put a kill button so user can close application easily 
st.write("""
**Once you're done hit the *Stop Server* button wait 5 seconds close the tab**
""")
if st.button("Stop Server"):
    st.warning("Stopping the Streamlit server...")
    os.kill(os.getpid(), signal.SIGTERM)
