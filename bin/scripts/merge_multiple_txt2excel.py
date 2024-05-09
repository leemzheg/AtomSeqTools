#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys
import argparse, textwrap
import pandas as pd
from natsort import natsorted

def argparse_line():
    parser = argparse.ArgumentParser(description='', 
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-txt', metavar='TXT', nargs='+', 
        help='multiple text files', required=True)
    parser.add_argument('-excel', metavar='EXCEL', 
        help='output excel file', required=True)
    argv = vars(parser.parse_args())
    return argv

def txt2excel(txt_files, excel):
    txt_list = []
    sheet_list = []
    for txt in natsorted(txt_files):
        sheet_list.append(txt.split('/')[-1].split('.')[0])
        txt_list.append(txt)

    writer = pd.ExcelWriter(excel, engine='xlsxwriter')
    pd.io.formats.excel.ExcelFormatter.header_style = None
    for i in range(0, len(txt_list)):
        df = pd.read_csv('{0}'.format(txt_list[i]), sep='\t')
        df.to_excel(writer, sheet_name='{0}'.format(sheet_list[i]), 
                    index=False)
        workbook = writer.book
        worksheet = writer.sheets[sheet_list[i]]
        header_format = workbook.add_format({
            'valign':'vcenter',
            'align':'center',
            'font_size':12,
            'font_name':'Arial',
            'border':1
            })
        default_format = workbook.add_format({
            'valign':'vcenter',
            'align':'center',
            'font_size':12,
            'font_name':'Arial',
            })
        none_center = workbook.add_format({
            'valign':'vcenter',
            'font_size':12,
            'font_name':'Arial',
            })
        worksheet.set_row(0, 30, header_format)
        worksheet.set_column('A:A', 30, default_format)
        worksheet.set_column('B:B', 13, default_format)
        worksheet.set_column('C:D', 18, default_format)
        worksheet.set_column('E:E', 16, default_format)
        worksheet.set_column('F:F', 10, default_format)
        worksheet.set_column('G:G', 18, default_format)
        worksheet.set_column('H:H', 18, none_center)
        worksheet.set_column('I:I', 10, default_format)
        worksheet.set_column('J:J', 18, default_format)
        worksheet.set_column('K:K', 18, none_center)
        worksheet.set_column('L:L', 20, none_center)
        worksheet.set_default_row(30)
    writer.save()

def main():
    argv = argparse_line()
    txt2excel(argv['txt'], argv['excel'])
if __name__ == '__main__':
    main()
