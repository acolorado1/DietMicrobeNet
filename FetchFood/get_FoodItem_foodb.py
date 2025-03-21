import streamlit as st 
import pandas as pd 
import os

# read in food information 
food_df = pd.read_csv('/Users/burkhang/Code_Projs/DietMicrobeNet/Data/food.csv')

# title of application
st.title("Get Metabolomes of Food Items Using FooDB")

# Search functionality
st.subheader("Search by Name")
search_query = st.text_input("Enter a name to search:")
if search_query:
    search_results = food_df[food_df["name"].str.contains(search_query, case=False)]
    st.write(search_results)
    if search_results.empty:
        st.warning("No matching results found!")

# Allow the user to select rows by index
selected_rows = st.multiselect(
    "Select rows by index",
    food_df.index.tolist(),  # Options are the index of the DataFrame
    default=[]  # Default selection is empty
)

# Show selected rows
if selected_rows:
    st.write("### Selected Rows")
    st.write(food_df.loc[selected_rows])  # Use .loc to pick rows by index
else:
    st.write("No rows selected.")

# put a kill button so user can close application easily 
st.write("""
**Once you're done hit the *Stop Server* button wait 5 seconds close the tab**
""")
if st.button("Stop Server"):
    st.warning("Stopping the Streamlit server...")
    os.kill(os.getpid(), signal.SIGTERM)