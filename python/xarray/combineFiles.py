# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 11:29:05 2024

@author: zhixi
"""

import os
from PyPDF2 import PdfMerger

def merge_pdfs(input_dir, output_filename='merged.pdf'):
    merger = PdfMerger()

    # Iterate through all PDF files in the input directory
    for i, filename in enumerate(os.listdir(input_dir)):
        if filename.endswith('.pdf'):
            filepath = os.path.join(input_dir, filename)
            merger.append(filepath)
            if i%50 == 0:
                print(f'merging file {i}')

    # Write the merged PDF to the output file
    with open(output_filename, 'wb') as output_file:
        merger.write(output_file)

    print(f"PDF files in '{input_dir}' merged into '{output_filename}' successfully.")
#%%

input_dir = r'C:\Users\zhixi\Downloads\allUnits\append'
output = f'{input_dir}\\merged.pdf'
merge_pdfs(input_dir, output_filename=output)
#%%

    