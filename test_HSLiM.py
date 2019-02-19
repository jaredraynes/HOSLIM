import pytest
import pandas as pd
import H_SLiM

def test_sequence_download1():
	download_output = H_SLiM.fasta_download("Q16082", 0, 65)

	expected = ('MSGRSVPHAHPATAEYEFANPSRLGEQRFGEGLLPEEILTPTLYHGYYVRPRAAPAGEGSRAGAS', 'Q16082', 'HSPB2', 'Human')

	assert download_output == expected

def test_HSLiM_output1():
	seq = H_SLiM.fasta_download('B6VPY2', -1, -1)
	scale = H_SLiM.scales("K_D")
	IDP = H_SLiM.order("IDP")
	HSLiM_output =  H_SLiM.h_slim(seq[0], seq[1], seq[2], seq[3], scale, IDP, 25, "", "no_save")

	data = [('as2-casein', 'B6VPY2', 'Domestic water buffalo', 10.3, 4, 3, 'FFIF', 6),
			('', '', '', 14.5, 8, 8, 'CLLAVALA', 15),
			('', '', '', 5.2, 3, 41, 'MAI', 43),
			('', '', '', 4.8, 3, 113, 'YLY', 115),	
			('', '', '', 8.2, 3, 119, 'IVL', 121)
			]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

def test_HSLiM_output2():
	seq = H_SLiM.fasta_download('D0QJ96', -1, -1)
	scale = H_SLiM.scales("K_D")
	IDP = H_SLiM.order("IDP")
	HSLiM_output =  H_SLiM.h_slim(seq[0], seq[1], seq[2], seq[3], scale, IDP, 25, "", "no_save")

	data = [('CSN1 ECO0000313EMBLACU257801', 'D0QJ96', 'Duckbill platypus', 29.0, 14, 3, 'VLILACLVAVAVAM', 16),
			('', '', '', 5.1, 3, 44, 'YYL', 46),
			('', '', '', 8.1, 3, 69, 'LLL', 71),
			('', '', '', 7.0, 3, 138, 'YFI', 140),
			('', '', '', 5.4, 4, 143, 'AAVY', 146),
			('', '', '', 6.6, 3, 151, 'LVY', 153),
			('', '', '', 4.5, 3, 169, 'YAF', 171)
			]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

def test_HSLiM_output3():
	seq = H_SLiM.fasta_download('D0QJA4', -1, -1)
	scale = H_SLiM.scales("K_D")
	IDP = H_SLiM.order("IDP")
	HSLiM_output =  H_SLiM.h_slim(seq[0], seq[1], seq[2], seq[3], scale, IDP, 25, "", "no_save")

	data = [('D0QJA4_9MAMM', 'D0QJA4', 'Australian echidna', 28.5, 14, 3, 'VFIFACLVAVAMAV', 16),
			('', '', '', 6.8, 3, 33, 'LVM', 35),
			('', '', '', 6.1, 3, 43, 'ALV', 45),
			('', '', '', 5.4, 3, 62, 'MVY', 64),
			('', '', '', 9.5, 4, 74, 'YIFF',77)
			]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

def test_HSLiM_output4():
	seq = H_SLiM.fasta_download('O97944', -1, -1)
	scale = H_SLiM.scales("K_D")
	IDP = H_SLiM.order("IDP")
	HSLiM_output =  H_SLiM.h_slim(seq[0], seq[1], seq[2], seq[3], scale, IDP, 25, "", "no_save")

	data = [('as2-casein', 'O97944', 'Dromedary', 9.9, 4, 3, 'FFIF', 6),
			('', '', '', 15.7, 8, 8, 'CLLAVVLA', 15),
			('', '', '', 6.0, 3, 42, 'VAI', 44),
			('', '', '', 6.7, 3, 98, 'IVM', 100),
			('', '', '', 6.5, 3, 167, 'FLW', 169)
			]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

def test_HSLiM_output5():
	seq = H_SLiM.fasta_download('Q9GKK3', -1, -1)
	scale = H_SLiM.scales("K_D")
	IDP = H_SLiM.order("IDP")
	HSLiM_output =  H_SLiM.h_slim(seq[0], seq[1], seq[2], seq[3], scale, IDP, 25, "", "no_save")

	data = [('beta-casein', 'Q9GKK3', 'Horse', 29.3, 13, 3, 'ILILACLVALALA', 15),
			('', '', '', 6.9, 3, 80, 'VVY', 82),
			('', '', '', 7.7, 4, 90, 'YAVV', 93),
			('', '', '', 7.2, 3, 176, 'LML', 178),
			('', '', '', 10.4, 5, 210, 'AFLLY', 214),
			('', '', '', 9.8, 4, 232, 'IVAV', 235),
			('', '', '', 9.1, 3, 239, 'VIV', 241)
			]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

def test_HSLiM_output6():
	seq = H_SLiM.fasta_download('B7VGF9', -1, -1)
	scale = H_SLiM.scales("K_D")
	IDP = H_SLiM.order("IDP")
	HSLiM_output =  H_SLiM.h_slim(seq[0], seq[1], seq[2], seq[3], scale, IDP, 25, "", "no_save")

	data = [('as2-casein', 'B7VGF9', 'Donkey', 10.4, 4, 3, 'FFIF', 6),
			('', '', '', 14.6, 8, 8, 'CLLAVALA', 15),
			('', '', '', 9.4,	4, 41, 'YVVI', 44),
			('', '', '', 8.2,	3, 130, 'IVL', 132)
			]
	
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

def test_HSLiM_output7():
	seq = H_SLiM.fasta_download('P02666', -1, -1)
	scale = H_SLiM.scales("K_D")
	IDP = H_SLiM.order("IDP")
	HSLiM_output =  H_SLiM.h_slim(seq[0], seq[1], seq[2], seq[3], scale, IDP, 25, "", "no_save")

	data = [('beta-casein', 'P02666', 'Bovine', 28.6, 13, 3, 'VLILACLVALALA', 15),
			('', '', '', 6.9, 3, 73, 'LVY', 75),
			('', '', '', 8.4, 3, 97, 'VVV', 99),
			('', '', '', 3.1, 3, 116, 'AMA', 118),
			('', '', '', 7, 3, 170, 'VMF', 172),
			('', '', '', 10.3, 5, 204, 'AFLLY', 208),
			('', '', '', 9.6, 3, 222, 'IIV', 224)
			]
	
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

	#ADD TEST FOR P06796

#THIS WORKS?? WTF beta-casein	Q9GKK3	Horse	29.3	13	3	ILILACLVALALA	15
