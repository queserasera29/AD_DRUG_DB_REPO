import streamlit as st
import pandas as pd

# Load data
@st.cache_data
def load_data():
    df = pd.read_csv("PubChem_compound_text_antidepressant drugs.csv")
    df.fillna("N/A", inplace=True)
    return df

df = load_data()

# Set page config
st.set_page_config(page_title="Antidepressant Drug Info", layout="wide", page_icon="💊")

# Custom CSS for aesthetics
st.markdown("""
    <style>
        body { font-family: 'Segoe UI', sans-serif; }
        .big-title { font-size: 2.5em; font-weight: bold; color: #00adb5; }
        .subtext { color: #aaaaaa; font-size: 1.1em; margin-bottom: 20px; }
        .stButton > button { background-color: #00adb5; color: white; }
    </style>
""", unsafe_allow_html=True)

# Landing Page
st.markdown('<div class="big-title">Antidepressant Drug Information Portal 💊</div>', unsafe_allow_html=True)
st.markdown('<div class="subtext">Explore chemical and pharmacological data on antidepressants using a simple search interface. Powered by PubChem data.</div>', unsafe_allow_html=True)

col1, col2 = st.columns([3, 1])
with col1:
    query = st.text_input("🔍 Enter Drug Name or Synonym")
with col2:
    drug_classes = sorted(df["Drug_Class"].dropna().unique())
    selected_class = st.selectbox("📂 Filter by Drug Class", options=["All"] + drug_classes)

if st.button("Search"):
    if query or selected_class != "All":
        filtered_df = df[
            (df["Name"].str.contains(query, case=False, na=False) | 
             df["Synonyms"].str.contains(query, case=False, na=False))
        ]
        if selected_class != "All":
            filtered_df = filtered_df[filtered_df["Drug_Class"] == selected_class]
        
        if filtered_df.empty:
            st.warning("No matching drug found.")
        else:
            for idx, row in filtered_df.iterrows():
                with st.expander(f"🧪 {row['Name']}"):
                    st.markdown(f"**Molecular Formula**: {row['Molecular_Formula']}")
                    st.markdown(f"**Molecular Weight**: {row['Molecular_Weight']}")
                    st.markdown(f"**XLogP**: {row['XLogP']}")
                    st.markdown(f"**Drug Class**: {row['Drug_Class']}")
                    
                    # Show PubChem image if CID exists
                    if row.get('Compound_CID') and row['Compound_CID'] != "N/A":
                        cid = str(row['Compound_CID']).strip()
                        img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/PNG"
                        st.image(img_url, caption="Structure from PubChem")

                    # Synonyms
                    st.markdown("**Synonyms:**")
                    if len(row["Synonyms"]) > 100:
                        st.markdown(row["Synonyms"][:100] + "... " + f"[View More](https://pubchem.ncbi.nlm.nih.gov/compound/{cid})")
                    else:
                        st.write(row["Synonyms"])
    else:
        st.warning("Please enter a drug name or select a drug class.")
