all:
	# [...
	python plot_XPS_C185_ref.py -r ../data/Data_XPS_C185_exp.csv \
	../data/results_C185/Data_XPS_C185_SJ_0.csv ../data/results_C185/Data_XPS_C185_SJ_Ef.csv ../data/results_C185/Data_XPS_C185_SJ_Ar.csv ../data/results_C185/Data_XPS_C185_SJ_Evac.csv ../data/results_C185/Data_XPS_C185_SJ_WF.csv \
	../data/results_C185/Data_XPS_C185_SJn_0.csv ../data/results_C185/Data_XPS_C185_SJn_Ef.csv ../data/results_C185/Data_XPS_C185_SJn_Ar.csv ../data/results_C185/Data_XPS_C185_SJn_Efp.csv \
	-n "SJ/\$$E_{ref}=0\$$" "SJ/\$$E_{ref}=E_F\$$" "SJ/\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ/\$$E_{ref}=E_\infty\$$" "SJ/\$$E_{ref}=\\phi\$$" \
	"SJn/\$$E_{ref}=0\$$" "SJn/\$$E_{ref}=E_F\$$" "SJn/\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJn/\$$E_{ref}=E_F'\$$" \
	-o ../Figure3.pdf
	# ...]
	# [...
	python plot_XPS_C185_2.py  -r ../data/Data_XPS_C185_exp.csv \
	../data/results_C185/Data_XPS_C185_SJn_0.csv ../data/results_C185/Data_XPS_C185_SJn_Ar.csv -n '$$E_{ref}=0$$' '$$E_{ref}=\varepsilon_{Ar,2s}$$' -o ../Figure4.pdf
	# ...]
	# SI
	python plot_surface_energy.py -o ../FigureS2.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C2H6 -n 'a (C$$_2$$H$$_6$$)' -o ../FigureS4a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C3H8 -n 'b (C$$_3$$H$$_8$$)' -o ../FigureS4b.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C2H4 -n 'c (C$$_2$$H$$_4$$)' -o ../FigureS5a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C3H6 -n 'd (C$$_3$$H$$_6$$)' -o ../FigureS5b.pdf
	python ./plot_XPS_slabs_adsorbate.py -a CO2 -n 'e (CO$$_2$$)' -o ../FigureS6a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a CH2O -n 'f (CH$$_2$$O)' -o ../FigureS6b.pdf
	python ./plot_XPS_slabs_adsorbate.py -a C2H4O -n 'g (C$$_2$$H$$_4$$O)' -o ../FigureS7a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a THF -n 'h (THF)' -o ../FigureS7b.pdf
	python ./plot_XPS_slabs_adsorbate.py -a CH4O -n 'i (CH$$_4$$O)' -o ../FigureS8a.pdf
	python ./plot_XPS_slabs_adsorbate.py -a THFn -n 'j (THF$$_n$$)' -o ../FigureS8b.pdf
	python ./plot_XPS_slabs_adsorbate_any.py -i ../data/Data_XPS_slab_adsorbate_SJn.csv -t SJn  -o ../FigureS9.pdf
clean:
	rm ../FigureS*.pdf ../Figure{3,4}.pdf
