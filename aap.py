# app.py
import streamlit as st
import pandas as pd
import base64
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO

# Load Data
df = pd.read_csv("PubChem_compound_text_antidepressant drugs.csv")
df.fillna("N/A", inplace=True)

# Extract Drug Classes from 'Tagged_by_PubChem'
def extract_class(tag):
    if tag == "N/A":
        return "Unknown"
    return tag.split(">")[-1].strip()

df["Drug_Class"] = df["Tagged_by_PubChem"].apply(extract_class)

# App Config
st.set_page_config(page_title="Antidepressant Drug Info Portal", layout="wide", initial_sidebar_state="auto")

# CSS Styling
dark_theme = """
<style>
body {
    background-color: #121212;
    color: #e0e0e0;
}
.stButton > button {
    border-radius: 8px;
    background: linear-gradient(to right, #6a11cb, #2575fc);
    color: white;
    transition: 0.3s;
}
.stButton > button:hover {
    transform: scale(1.05);
    box-shadow: 0px 0px 10px rgba(255,255,255,0.2);
}
</style>
"""
st.markdown(dark_theme, unsafe_allow_html=True)

# Title and Intro
st.title("\U0001F48A Antidepressant Drug Information Portal")
st.markdown("""
Welcome to the Antidepressant Drug Information Portal — a sleek, searchable tool for exploring chemical and pharmacological details of antidepressants.
Use the search bar below to find compounds by name or synonym, or filter by drug class.
""")

# Search Interface
col1, col2, col3 = st.columns([3, 2, 1])
with col1:
    query = st.text_input("Enter Drug Name or Synonym:")
with col2:
    drug_classes = ["All"] + sorted(df["Drug_Class"].unique())
    selected_class = st.selectbox("Filter by Drug Class:", drug_classes)
with col3:
    search_triggered = st.button("\U0001F50D Search")

# Filter Logic
results = pd.DataFrame()
if search_triggered:
    q = query.lower().strip()
    results = df[df["Name"].str.lower().str.contains(q) | df["Synonyms"].str.lower().str.contains(q)]
    if selected_class != "All":
        results = results[results["Drug_Class"] == selected_class]

    if not results.empty:
        st.subheader("\U0001F4C4 Search Results")
        for idx, row in results.iterrows():
            with st.expander(f"{row['Name']} ({row['Molecular_Formula']})"):
                st.markdown(f"**Molecular Weight**: {row['Molecular_Weight']}")
                st.markdown(f"**XLogP**: {row['XLogP']}")
                st.markdown(f"**Drug Class**: {row['Drug_Class']}")

                short_synonyms = row['Synonyms'].split('|')[:3]
                full_synonyms = row['Synonyms'].replace('|', ', ')
                st.markdown(f"**Synonyms**: {', '.join(short_synonyms)} ...")
                with st.expander("View all synonyms"):
                    st.code(full_synonyms)

                # Draw molecule
                if row['SMILES'] != "N/A":
                    mol = Chem.MolFromSmiles(row['SMILES'])
                    if mol:
                        img = Draw.MolToImage(mol, size=(300, 300))
                        st.image(img, caption="Structure")

                # Download option
                def convert_df_to_csv(df):
                    return df.to_csv(index=False).encode('utf-8')

                csv_data = convert_df_to_csv(pd.DataFrame([row]))
                st.download_button("Download Drug Info (CSV)", csv_data, file_name=f"{row['Name']}.csv", mime='text/csv')
    else:
        st.error("No matching drugs found. Please try another name or class.")

# Footer
st.markdown("""
---
App built using **Streamlit**, **RDKit**, and **Pandas**.
Inspired by **PubChem** but tailored for focused antidepressant data.
""")
