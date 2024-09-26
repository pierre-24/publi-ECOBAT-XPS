all:
	python plot_XPS_C185_2.py -i ../data/Data_XPS_C185_SJn_0.csv -o ../Figure4.pdf
	python plot_XPS_trends.py -o ../Figure5.pdf
	python plot_XPS_slabs.py -o ../Figure6.pdf
	python plot_XPS_slabs_adsorbate_any.py -o ../Figure12.pdf
	python plot_surface_energy.py -o ../FigureA1.pdf
	python plot_XPS_C185.py -i ../data/Data_XPS_C185_fermi.csv -o ../FigureA2.pdf
	python plot_XPS_slabs.py -i ../data/Data_XPS_slab_SJn.csv -t SJn -o ../FigureA3.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C2H6 -n 'a (C$$_2$$H$$_6$$)' -o ../FigureA4a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C3H8 -n 'b (C$$_3$$H$$_8$$)' -o ../FigureA4b.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C2H4 -n 'c (C$$_2$$H$$_4$$)' -o ../FigureA5a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C3H6 -n 'd (C$$_3$$H$$_6$$)' -o ../FigureA5b.pdf
	python ./plot_XPS_slabs_adsorbate.py -a CO2 -n 'e (CO$$_2$$)' -o ../FigureA6a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a CH2O -n 'f (CH$$_2$$O)' -o ../FigureA6b.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C2H4O -n 'g (C$$_2$$H$$_4$$O)' -o ../FigureA7a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a THF -n 'h (THF)' -o ../FigureA7b.pdf
	python ./plot_XPS_slabs_adsorbate.py -a CH4O -n 'i (CH$$_4$$O)' -o ../FigureA8a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a THFn -n 'j (THF$$_n$$)' -o ../FigureA8b.pdf
	python ./plot_XPS_slabs_adsorbate_any.py -i ../data/Data_XPS_slab_adsorbate_SJn.csv -t SJn  -o ../FigureA9.pdf
	# tmp
	python plot_XPS_C185_ref.py -r ../data/Data_XPS_C185_exp.csv ../data/Data_XPS_C185_SJ_{0,Ef,Ar,WF}.csv ../data/Data_XPS_C185_SJn_{0,Ef,Ar}.csv -n "SJ/\$E_{ref}=0\$" "SJ/\$E_{ref}=E_F\$" "SJ/\$E_{ref}=\\varepsilon_{Ar}\$" "SJ/\$E_{ref}=\$WF" "SJn/\$E_{ref}=0\$" "SJn/\$E_{ref}=E_F\$" "SJn/\$E_{ref}=\\varepsilon_{Ar}\$"  -o ../Figure4A.pdf
