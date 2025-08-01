import streamlit as st
import pandas as pd

@st.cache_data
def load_data():
    return pd.read_csv("PubChem_compound_text_antidepressant drugs.csv").fillna("N/A")

df = load_data()

st.set_page_config(page_title="Antidepressant Drug Info", layout="wide")
st.title("Antidepressant Drug Information Portal")

query = st.text_input("Enter Drug Name or Synonym")
selected = st.selectbox("Filter by Drug Class", options=["All"] + sorted(df["Name"].unique()))

if st.button("Search"):
    filtered = df[df["Name"].str.contains(query, case=False, na=False) |
                  df["Synonyms"].str.contains(query, case=False, na=False)]
    if selected != "All":
        filtered = filtered[filtered["Drug_Class"] == selected]

    if filtered.empty:
        st.warning("No results found")
    else:
        st.dataframe(filtered[["Name", "Molecular_Formula", "Molecular_Weight", "XLogP"]])
