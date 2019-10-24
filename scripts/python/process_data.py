import pandas as pd

BASES = ['A', 'C', 'G', 'T']
GAP = '-'

def get_state(wildtype, genotype):
	if genotype == GAP:
		return '?'
	elif genotype == wildtype:
		return '0'
	elif genotype not in BASES:
		return '1'
	else:
		return '2'

def trim_df(df):
	df.columns = df.columns.str.strip()
	for col in df.columns:
		df[col] = df[col].str.strip()
	return df

# transform mutation data to ternary format
wildtype = 'LN-T1'
extra_cols = ['Coordinates', 'LN-T1', 'LC-T1']
filename = 'mutations.csv'

df = trim_df(pd.read_csv(filename))
data_cols = df.columns.difference(extra_cols)
for col in data_cols:
	df[col] = df.apply(lambda row: get_state(row[wildtype], row[col]), axis=1)

res = df.agg(lambda x: ''.join(x))
res_df = pd.DataFrame(res, res.index, ['sequence'])
res_df.drop(extra_cols, inplace=True)
res_df.to_csv('single_cell.csv', index_label='cell')

# create beast alignment xml objects
res_df = pd.read_csv('single_cell.csv')

align_start = "<data id='alignment' dataType='ternaryWithError'>\n"
align_end = "\n</data>"
seq_start = "\t<sequence taxon='%s'>\n\t %s"
seq_end = "\n\t</sequence>"

res_df['xml'] = res_df.apply(lambda x: seq_start % (x['cell'], x['sequence']) + seq_end, axis=1)
xml_lines = res_df['xml'].apply(lambda x: ''.join(x))
xml_lines = align_start + '\n'.join(xml_lines) + align_end

outfile = open('single_cell.xml', 'w')
outfile.writelines(xml_lines)
outfile.close()