import PyPDF2
from PIL import Image
import io

def image_to_pdf_page(image_path):
    """Converts an image file into a PyPDF2 page object."""
    img = Image.open(image_path)
    
    # Convert RGBA to RGB (PDF doesn't support alpha channels directly this way)
    if img.mode == 'RGBA':
        img = img.convert('RGB')
        
    pdf_bytes = io.BytesIO()
    img.save(pdf_bytes, format='PDF', resolution=100.0)
    pdf_bytes.seek(0)
    
    return PyPDF2.PdfReader(pdf_bytes).pages[0]

def assemble_hw_pdf(original_pdf, output_pdf, insertions):
    """Iterates through the PDF and inserts image pages at specified indices."""
    print("Reading original PDF...")
    reader = PyPDF2.PdfReader(original_pdf)
    writer = PyPDF2.PdfWriter()

    # Iterate through all pages in the handwritten PDF
    for i in range(len(reader.pages)):
        # 1. Add the handwritten page
        writer.add_page(reader.pages[i])
        
        # 2. Check if we need to insert plots AFTER this page
        # Note: PyPDF2 is 0-indexed, so Page 1 is index 0.
        if i in insertions:
            for img_path in insertions[i]:
                print(f"Inserting {img_path} after Page {i + 1}...")
                plot_page = image_to_pdf_page(img_path)
                writer.add_page(plot_page)

    # Save the final compiled document
    print(f"Saving compiled document to {output_pdf}...")
    with open(output_pdf, "wb") as f:
        writer.write(f)
    print("Done!")

# --- Execution ---
# Dictionary format: { page_index_to_insert_after: ["list", "of", "images"] }
# Remember: Index is Page Number - 1
plot_insertions = {
    1: ["levenspiel_plot.png", "pfr_volume_vs_conversion.png"], # After Page 2
    3: ["multistage analysis.png"],                             # After Page 4
    7: ["coolant.png", "coolant_b.png"],                        # After Page 8
    11: ["g(T)_R(T).png", "g(T)_R(T)_no_coolant.png"]           # After Page 12
}

assemble_hw_pdf(
    original_pdf="hw5 (1).pdf", 
    output_pdf="hw5_final_submission.pdf", 
    insertions=plot_insertions
)