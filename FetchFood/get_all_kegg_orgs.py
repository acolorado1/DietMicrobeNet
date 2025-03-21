import streamlit as st
import requests

# this code was written with the help of chatGPT 
def fetch_kegg_organisms():
    """
    Fetch the full list of organisms from KEGG.
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
            print(parts)
            if len(parts) == 4:
                t_number, code, scientific_name, taxonomy = parts
                organisms.append((code, scientific_name, taxonomy))
        return organisms
    else:
        raise Exception("Failed to fetch data from KEGG.")

def search_organism_by_name(query, organisms):
    """
    Search for organisms by name (scientific or common).
    Args:
        query (str): Search term.
        organisms (list): List of organisms fetched from KEGG.
    Returns:
        List of matching organisms.
    """
    query = query.lower()
    results = [
        org for org in organisms
        if query in org[1].lower() or query in org[2].lower()
    ]
    return results

# Streamlit App
st.title("KEGG Organism Search")

st.write("""
This app allows you to search for organisms in the KEGG database by their scientific or common names. 
You can also view the KEGG codes associated with these organisms.
""")

# Load KEGG organism data
with st.spinner("Fetching KEGG organism data..."):
    try:
        organisms = fetch_kegg_organisms()
        st.success(f"Successfully loaded {len(organisms)} organisms!")
    except Exception as e:
        st.error("Failed to fetch KEGG data. Please try again later.")
        st.stop()

# Search bar
query = st.text_input("Enter an organism name to search (scientific or common):", "")

# Results display
if query:
    matches = search_organism_by_name(query, organisms)
    if matches:
        st.write(f"Found {len(matches)} result(s):")
        for match in matches:
            st.write(f"- **Code**: `{match[0]}` | **Scientific Name**: *{match[1]}* | **Common Name**: {match[2]}")
    else:
        st.write("No matching organisms found. Please try again.")

# Optional: Display all KEGG codes
if st.checkbox("Show all KEGG organism codes"):
    codes = [org[0] for org in organisms]
    st.write(", ".join(codes))
