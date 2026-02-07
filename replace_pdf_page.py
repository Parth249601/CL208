from PyPDF2 import PdfReader, PdfWriter
import os

# Define file paths
main_pdf_path = "hw2/CRE-HW2-Part1.pdf"
replacement_pdf_path = "hw2/q2.pdf"
output_pdf_path = "hw2/CRE-HW2-Part1_updated.pdf"

# Read the main PDF
main_reader = PdfReader(main_pdf_path)
output_writer = PdfWriter()

# Read the replacement PDF
replacement_reader = PdfReader(replacement_pdf_path)

# Get the first page from q2.pdf to replace page 3 (index 2)
replacement_page = replacement_reader.pages[0]

# Copy all pages from the main PDF, replacing page 3
for page_num in range(len(main_reader.pages)):
    if page_num == 2:  # Page 3 is at index 2 (0-indexed)
        output_writer.add_page(replacement_page)
    else:
        output_writer.add_page(main_reader.pages[page_num])

# Write the output PDF
with open(output_pdf_path, 'wb') as output_file:
    output_writer.write(output_file)

print(f"Successfully replaced page 3!")
print(f"Output saved to: {output_pdf_path}")

# Optional: Replace the original file
import shutil
shutil.copy(output_pdf_path, main_pdf_path)
print(f"Original file updated: {main_pdf_path}")
