# ...
function waread3_py(x, flag=true)
    #     # use py"""...""" syntax to access the Python function
    #     py"""
    #     import pandas as pd
    #     import datetime

    #     def waread3(x, flag=True):
    #         if flag:
    #             df = pd.read_csv(x, delim_whitespace=True, header=0,
    #                             na_values=-9999, verbose=True,engine='c')
    #             if 'YY' not in df.columns:
    #                 print("Column 'YY' not found in the CSV file.")
    #                 return None
    #             if df.iloc[:, 0].dtype == 'object':
    #                 print('First column contains strings, subsetting to Int...')
    #                 df = df[~df.iloc[:, 0].str.contains("[A-z]|-", na=False)]
    #             source_col_loc = df.columns.get_loc('YY')        
    #             df['date'] = df.iloc[:, source_col_loc:source_col_loc +
    #                                 3].apply(lambda x: "-".join(x.astype(str)), axis=1)
    #             df = df.iloc[:, 4:]
    #             df['date'] = pd.to_datetime(df['date'])
    #             #df.set_index('date', inplace=True)
    #             df.iloc[:,0:-2] = df.iloc[:,0:-2].apply(lambda x: x.astype(float), axis=1)
    #             df.filename = x
    #             print(df.filename,"done!")
    #             return df
    #         else:
    #             print('Date parse failed, try reading without date transformation...')
    #             df = pd.read_csv(x, delim_whitespace=True, comment="Y", skip_blank_lines=True).dropna()
    #             df.filename = x
    #             print(df.filename,"done!")
    #             return df
    #     """
    #     # call the Python 