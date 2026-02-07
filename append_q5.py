from PyPDF2 import PdfReader, PdfWriter

# Define file paths
main_pdf_path = "hw2/CRE-HW2-Part1_updated.pdf"
append_pdf_path = "hw2/q5.pdf"
output_pdf_path = "hw2/CRE-HW2-Part1_updated.pdf"

# Read the main PDF
main_reader = PdfReader(main_pdf_path)
output_writer = PdfWriter()

# Read the PDF to append
append_reader = PdfReader(append_pdf_path)

# Add all pages from the main PDF
for page in main_reader.pages:
    output_writer.add_page(page)

# Add all pages from q5.pdf
for page in append_reader.pages:
    output_writer.add_page(page)

# Write the output PDF (overwriting the updated file)
with open(output_pdf_path, 'wb') as output_file:
    output_writer.write(output_file)

print(f"Successfully appended q5.pdf!")
print(f"Updated file: {output_pdf_path}")
print(f"Total pages: {len(output_writer.pages)}")
