all:
	# [...
	python plot_XPS_C185_ref.py -r ../data/Data_XPS_C185_exp.csv \
	../data/results_C185/Data_XPS_C185_SJ_0.csv ../data/results_C185/Data_XPS_C185_SJ_Ef.csv ../data/results_C185/Data_XPS_C185_SJ_Ar.csv ../data/results_C185/Data_XPS_C185_SJ_Evac.csv ../data/results_C185/Data_XPS_C185_SJ_WF.csv \
	../data/results_C185/Data_XPS_C185_SJn_0.csv ../data/results_C185/Data_XPS_C185_SJn_Ef.csv ../data/results_C185/Data_XPS_C185_SJn_Ar.csv ../data/results_C185/Data_XPS_C185_SJn_Efp.csv \
	-n "SJ/\$$E_{ref}=0\$$" "SJ/\$$E_{ref}=E_F\$$" "SJ/\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ/\$$E_{ref}=E_\infty\$$" "SJ/\$$E_{ref}=\\phi\$$" \
	"SJ\$$^n$$/\$$E_{ref}=0\$$" "SJ\$$^n$$/\$$E_{ref}=E_F\$$" "SJ\$$^n$$/\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ\$$^n$$/\$$E_{ref}=E_F'\$$" \
	-o ../Figure3.pdf
	# ...]
	# [...
	python plot_XPS_C185_2.py  -r ../data/Data_XPS_C185_exp.csv \
	../data/results_C185/Data_XPS_C185_SJn_0.csv ../data/results_C185/Data_XPS_C185_SJn_Ar.csv -n 'SJ$$^n$$+$$E_{ref}=0$$' 'SJ$$^n$$+$$E_{ref}=\varepsilon_{Ar,2s}$$' -o ../Figure4.pdf
	# ...]
	python ./plot_XPS_slabs_ref_2.py ../data/results_slabs/Data_XPS_slabs_SJ_0.csv ../data/results_slabs/Data_XPS_slabs_SJ_Ef.csv ../data/results_slabs/Data_XPS_slabs_SJ_Ar.csv ../data/results_slabs/Data_XPS_slabs_SJ_Evac.csv ../data/results_slabs/Data_XPS_slabs_SJ_WF.csv \
	-n "SJ+\$$E_{ref}=0\$$" "SJ+\$$E_{ref}=E_F\$$" "SJ+\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ+\$$E_{ref}=E_\infty\$$" "SJ+\$$E_{ref}=\\phi\$$" \
	-o ../Figure5.pdf
	python ./plot_XPS_slabs_ref_2.py ../data/results_slabs/Data_XPS_slabs_SJn_0.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ef.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ar.csv ../data/results_slabs/Data_XPS_slabs_SJn_Efp.csv \
	-n "SJ\$$^n$$+\$$E_{ref}=0\$$" "SJ\$$^n$$+\$$E_{ref}=E_F\$$" "SJ\$$^n$$+\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ\$$^n$$+\$$E_{ref}=E_F'\$$" \
	-o ../Figure6.pdf
	python ./plot_XPS_slabs_2.py ../data/results_slabs/Data_XPS_slabs_SJ_0.csv ../data/results_slabs/Data_XPS_slabs_SJ_Ef.csv ../data/results_slabs/Data_XPS_slabs_SJ_Ar.csv ../data/results_slabs/Data_XPS_slabs_SJ_Evac.csv ../data/results_slabs/Data_XPS_slabs_SJ_WF.csv \
	-n "SJ+\$$E_{ref}=0\$$" "SJ+\$$E_{ref}=E_F\$$" "SJ+\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ+\$$E_{ref}=E_\infty\$$" "SJ+\$$E_{ref}=\\phi\$$" \
	 -a Ca:-3:5 O:-10:2 -o ../Figure7a.pdf
	python ./plot_XPS_slabs_2.py ../data/results_slabs/Data_XPS_slabs_SJn_0.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ef.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ar.csv ../data/results_slabs/Data_XPS_slabs_SJn_Efp.csv \
	-n "SJ\$$^n$$+\$$E_{ref}=0\$$" "SJ\$$^n$$+\$$E_{ref}=E_F\$$" "SJ\$$^n$$+\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ\$$^n$$+\$$E_{ref}=E_F'\$$" \
	 -a Ca:-3:5 O:-10:4 -o ../Figure7b.pdf
	# SI
	python plot_surface_energy.py -o ../FigureS2.pdf
	python ./plot_XPS_slabs_ref.py ../data/results_slabs/Data_XPS_slabs_SJ_0.csv ../data/results_slabs/Data_XPS_slabs_SJ_Ef.csv ../data/results_slabs/Data_XPS_slabs_SJ_Ar.csv ../data/results_slabs/Data_XPS_slabs_SJ_Evac.csv ../data/results_slabs/Data_XPS_slabs_SJ_WF.csv \
	-n "SJ+\$$E_{ref}=0\$$" "SJ+\$$E_{ref}=E_F\$$" "SJ+\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ+\$$E_{ref}=E_\infty\$$" "SJ+\$$E_{ref}=\\phi\$$" \
	-o ../FigureS3.pdf \
	-t TableS3.tab
	python ./plot_XPS_slabs_ref.py ../data/results_slabs/Data_XPS_slabs_SJn_0.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ef.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ar.csv ../data/results_slabs/Data_XPS_slabs_SJn_Efp.csv \
	-n "SJ\$$^n\$$+\$$E_{ref}=0\$$" "SJ\$$^n\$$+\$$E_{ref}=E_F\$$" "SJ\$$^n\$$+\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ\$$^n\$$+\$$E_{ref}=E_F'\$$" \
	-o ../FigureS4.pdf \
	-t TableS4.tab
	python ./plot_XPS_slabs_adsorbate_2.py -a C2H6 -n 'a (C$$_2$$H$$_6$$)' -o ../FigureS5a.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a C3H8 -n 'b (C$$_3$$H$$_8$$)' -o ../FigureS5b.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a C2H4 -n 'c (C$$_2$$H$$_4$$)' -o ../FigureS6a.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a C3H6 -n 'd (C$$_3$$H$$_6$$)' -o ../FigureS6b.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a CO2 -n 'e (CO$$_2$$)' -o ../FigureS7a.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a CH2O -n 'f (CH$$_2$$O)' -o ../FigureS7b.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a C2H4O -n 'g (C$$_2$$H$$_4$$O)' -o ../FigureS8a.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a THF -n 'h (THF)' -o ../FigureS8b.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a CH4O -n 'i (CH$$_4$$O)' -o ../FigureS9a.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a THFn -n 'j (THF$$_n$$)' -o ../FigureS9b.pdf -x "$$(cat slabs_annotate.txt | tr '\n'  ',')"
clean:
	rm ../FigureS*.pdf ../Figure{3,4,5,6}.pdf
