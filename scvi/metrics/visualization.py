import seaborn as sns

color_dictionary = dict()

macosko_regev_labels = [['Pvalb low', 'Pvalb', 'Pvalb 1', 'Pvalb 2'],
                        ['Pvalb Ex_1', 'Pvalb Ex_2', 'Pvalb Ex'],
                        ['Pvalb Astro_1', 'Pvalb Astro_2'],
                        ['L2/3 IT Astro', 'L2/3 IT Macc1', 'L2/3 IT Sla_Astro', 'L2/3 IT', 'L2/3 IT Sla',
                         'L2/3 IT Sla_Inh'],
                        ['Sst Tac2', 'Sst Myh8', 'Sst Etv1', 'Sst Chodl', 'Sst'],
                        ['L5 PT_2', 'L5 PT IT', 'L5 PT_1'],
                        ['L5 IT Tcap_1_3', 'L5 IT Tcap_2', 'L5 IT Tcap_Astro', 'L5 IT Tcap_1', 'L5 IT Tcap_L2/3',
                         'L5 IT Tcap_Foxp2', 'L5 IT Tcap_3'],
                        ['L5 IT Aldh1a7_2', 'L5 IT Aldh1a7', 'L5 IT Aldh1a7_1'],
                        ['L5 NP', 'L5 NP Slc17a8'],
                        ['L6 IT Car3', 'L6 CT Olig', 'L6 IT Maf', 'L6 IT Ntn5 Mgp', 'L6 IT Ntn5 Inpp4b'],
                        ['L6 CT Nxph2', 'L6 CT Astro', 'L6 CT', 'L6 CT Grp'],
                        ['L6b', 'L6b F2r'],
                        ['Lamp5 Sncg', 'Lamp5 Egln3', 'Lamp5 Slc35d3'],
                        ['Vip Rspo4', 'Vip Serpinf1', 'Vip'],
                        ['Astro Ex', 'Astro Aqp4'],
                        ['OPC Pdgfra'],
                        ['VLMC Osr1'],
                        ['Oligo Enpp6_1', 'Oligo Enpp6_2', 'Oligo Opalin'],
                        ['Sncg Ptprk'],
                        ['Endo Slc38a5', 'Endo Slc38a5_Peri_2', 'Endo Slc38a5_Peri_1']]

macosko_regev_colors = [sns.color_palette("Greens")[2:6],  # Pvalb
                        sns.light_palette("green", 5)[0:3],  # Pvalb Ex
                        sns.light_palette("green", 5)[3:5],  # Pvalb Astro
                        sns.light_palette("orange", 6),  # L2/3
                        sns.light_palette('red')[1:6],  # Sst
                        sns.light_palette("cyan", 3),  # L5 PT
                        sns.light_palette('purple', 8)[1:8],  # L5 IT Tcap
                        sns.light_palette('purple', 7)[4:7],  # L5 IT Aldh1a7
                        sns.light_palette("navy", 7)[3:5],  # L5 NP
                        sns.light_palette("brown", 7)[2:7],  # L6 IT
                        sns.dark_palette("brown", 8)[1:5],  # L6 CT
                        sns.dark_palette("green", 8)[5:7],  # L6
                        sns.dark_palette("yellow", 7)[1:4],  # Lamp5
                        sns.dark_palette("yellow", 7)[4:7],  # Vip
                        sns.color_palette("Paired", 4),  # Astro OPC VLMC
                        sns.color_palette('Greys', 3),  # Oligo
                        sns.dark_palette('tan'),  # sncg
                        sns.light_palette('hotpink', 3)]  # endo]

macosko_regev_labels_flattened = [l1 for l2 in macosko_regev_labels for l1 in l2]
macosko_regev_colors = [c1 for c2 in macosko_regev_colors for c1 in c2]

macosko_regev_dictionary = dict(zip(macosko_regev_labels_flattened, macosko_regev_colors))

color_dictionary.update(macosko_regev_dictionary)

pancreas_dictionary = {
    'alpha': sns.color_palette("Blues")[0],
    'beta': sns.color_palette("Greens")[0],
    'gamma': sns.color_palette("Reds")[0],
    'delta': sns.dark_palette("yellow", 3)[2],
    'co-expression': sns.color_palette("Greys")[0],
}

color_dictionary.update(pancreas_dictionary)
