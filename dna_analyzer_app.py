import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from collections import Counter
import matplotlib.pyplot as plt

st.title("Unique DNA Code Analyzer")

# ---- Sequence Input ----
st.header("Enter DNA Sequence")
sequence = st.text_area("Paste your DNA sequence (A, T, C, G only):").strip().upper()

if sequence:
    # Basic Quality Check
    if not all(base in "ATCG" for base in sequence):
        st.error("Please enter a valid DNA sequence (A, T, C, G only).")
    else:
        st.subheader("Analysis Modules")

        # ---- 1. GC Content ----
        gc = gc_fraction(sequence) * 100
        st.markdown(f"**GC Content:** {gc:.2f}%")

        # GC Content Visualization
        fig, ax = plt.subplots()
        ax.bar(["G", "C"], [sequence.count("G"), sequence.count("C")])
        ax.set_ylabel("Count")
        st.pyplot(fig)

        # ---- 2. Motif Search ----
        motif = st.text_input("Motif to search (e.g., ATG, TATA):").strip().upper()
        if motif:
            count = sequence.count(motif)
            positions = [i+1 for i in range(len(sequence)) if sequence.startswith(motif, i)]
            st.info(f"Found the motif '{motif}' {count} times at positions: {positions}")

        # ---- 3. Translation ----
        st.subheader("Protein Translation (first frame)")
        protein = Seq(sequence).translate(to_stop=True)
        st.code(str(protein))

        # ---- 4. Automated Report ----
        st.subheader("Summary Report")
        report = f"""
        Your sequence is {len(sequence)} base pairs long.
        The GC content is {gc:.2f}%, which {'suggests stability' if gc > 50 else 'is typical for many genes'}.
        {f"The motif '{motif}' appears {count} times." if motif else ""}
        The first predicted protein fragment (from DNA translation) is:
        {str(protein)}
        """
        st.text_area("Automated Report", value=report.strip(), height=200)

        # ---- Educational Mode ----
        if st.checkbox("Educational Mode: Show explanations"):
            st.markdown("- **GC Content:** Higher GC can mean more stable DNA.")
            st.markdown("- **Motif Analysis:** Motifs are patterns with potential biological significance.")
            st.markdown("- **Translation:** Converts DNA to protein (amino acid sequence).")
