from pdfrw import PdfReader, PdfWriter

def invert_pdf_colors(filepath):
  # Read the PDF file
  reader = PdfReader(filepath)

  # Invert the colors of each page in the PDF
  for page in reader.pages:
      if 'Contents' in page:
          for instruction in page.Contents.stream:
              # Invert the colors by negating the values of the RGB color space
              if instruction.operator == 'rg':
                  instruction.args = [-x for x in instruction.args]
              elif instruction.operator == 'RG':
                  instruction.args = [-x for x in instruction.args]

  # Write the inverted PDF file
  writer = PdfWriter()
  writer.write(filepath + '_inverted.pdf', reader)

arg ="book.pdf"

invert_pdf_colors(arg)