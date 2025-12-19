import os
import tabulate
import numpy as np
import glob

def get_errors(NJOY_outdir:str, nuclides:list[str]=['U235']):
    '''
    Docstring for get_errors
    
    :param NJOY_outdir: Description
    :type NJOY_outdir: str
    :param nuclides: Description
    :type nuclides: list[str] or 'all'
        If 'all', then find all files with .out extension


    '''


    if nuclides == 'all':
        nuclides = sorted(glob.glob(f'{NJOY_outdir}/*.out'))
        nuclides = [i.strip('/')[-1][0:-4] for i in nuclides]

    error_descrip = []

    for nuc in nuclides:

        assert os.path.exists( f'{NJOY_outdir}/{nuc}.out'), f"Requested nuclide {nuc} not found."

        print(f'Isotope: {nuc}')
        with open(f'{NJOY_outdir}/{nuc}.out','r') as f:
            lines = f.readlines()

        error_found = False
        for line in lines:
            if '***error' in line:
                error_found = True
                error_descrip.append(line)

        if error_found == False:
            error_descrip.append('ok')

    return nuclides, error_descrip

def write_only_errors_to_file(NJOY_outpath, nuclides:list[str]=['U235']):

    nuc_njoy_outputs, error_descrip = get_errors(NJOY_outpath, nuclides)

    content = np.zeros((len(nuc_njoy_outputs),2), dtype=object)

    content[:,0] = nuc_njoy_outputs
    content[:,1] = error_descrip

    tab = Table(caption='NJOY error descriptions', table_ref='error_descrip', headers=['Isotope','NJOY error'], content=content)

    table_out = tab.write_to_latex()
    
    return table_out

class Table:
    def __init__(self, caption:str, table_ref:str, headers:list, content:np.ndarray):
        self.headers = headers
        self.caption = caption
        self.table_ref = table_ref
        assert len(self.headers) == content.shape[1]
        self.content = content

    def write_to_latex(self) -> str:
        '''
        outputs string 

        Arguments
        ---------
        precision_float : if entry is a float, change precision and output to 1 \times 10^{3} format. If None, then ignore

        write using the following format:
        \begin{table}[H]
            \caption{Basic information regarding dataset}
            \begin{center}
                \begin{tabular}{<ncol>*c}
                \toprule
                <headers>
                \midrule
                <content>
                \bottomrule
                \end{tabular}
            \label{tab:VrecMARENnuclides}
            \end{center}
        \end{table}

        note modify from the output of tabulate


        '''

        content = self.content
        tmp = tabulate.tabulate(tabular_data=content, headers=self.headers, tablefmt='latex_raw')
        tmp = tmp.splitlines()
        
        header_line = [tmp[2]]
        content_lines = tmp[4:4+self.content.shape[0]]
        ncol = 'c'*len(self.headers)
        front_mat = ['\\begin{table}[H]', '\\caption{'+f'{self.caption}'+'}', '\\begin{center}', '\\begin{tabular}{' + ncol + '}', '\\toprule']

        final = front_mat + header_line + ['\\midrule']  + content_lines + ['\\bottomrule', '\\end{tabular}', '\\label{tab:' + f'{self.table_ref}' + '}' ,'\\end{center}','\\end{table}']
        final = "\n".join(final)

        return final
