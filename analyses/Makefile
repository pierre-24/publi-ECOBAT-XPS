all: main SI

main:
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
	 -a Ca:-2:4 O:-10:4  -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ',')" -o ../Figure7a.pdf
	python ./plot_XPS_slabs_2.py ../data/results_slabs/Data_XPS_slabs_SJn_0.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ef.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ar.csv ../data/results_slabs/Data_XPS_slabs_SJn_Efp.csv \
	-n "SJ\$$^n$$+\$$E_{ref}=0\$$" "SJ\$$^n$$+\$$E_{ref}=E_F\$$" "SJ\$$^n$$+\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ\$$^n$$+\$$E_{ref}=E_F'\$$" \
	 -a Ca:-2:4 O:-10:4 -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ',')" -o ../Figure7b.pdf

SI:
	python plot_surface_energy.py -o ../FigureS2.pdf
	python ./plot_XPS_slabs_ref.py ../data/results_slabs/Data_XPS_slabs_SJ_0.csv ../data/results_slabs/Data_XPS_slabs_SJ_Ef.csv ../data/results_slabs/Data_XPS_slabs_SJ_Ar.csv ../data/results_slabs/Data_XPS_slabs_SJ_Evac.csv ../data/results_slabs/Data_XPS_slabs_SJ_WF.csv \
	-n "SJ+\$$E_{ref}=0\$$" "SJ+\$$E_{ref}=E_F\$$" "SJ+\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ+\$$E_{ref}=E_\infty\$$" "SJ+\$$E_{ref}=\\phi\$$" \
	-o ../FigureS3.pdf \
	-t TableS3.tab
	python ./plot_XPS_slabs_ref.py ../data/results_slabs/Data_XPS_slabs_SJn_0.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ef.csv ../data/results_slabs/Data_XPS_slabs_SJn_Ar.csv ../data/results_slabs/Data_XPS_slabs_SJn_Efp.csv \
	-n "SJ\$$^n\$$+\$$E_{ref}=0\$$" "SJ\$$^n\$$+\$$E_{ref}=E_F\$$" "SJ\$$^n\$$+\$$E_{ref}=\\varepsilon_{Ar,2s}\$$" "SJ\$$^n\$$+\$$E_{ref}=E_F'\$$" \
	-o ../FigureS4.pdf \
	-t TableS4.tab
	python ./plot_XPS_slabs_adsorbate_2.py -a C2H6 -n 'a (H$$_3$$C-CH$$_3$$)' -o ../FigureS5a.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/C2H6_ads.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a C3H8 -n 'b (H$$_3$$C-CH$$_2$$-CH$$_3$$)' -o ../FigureS5b.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/C3H8_ads.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a C2H4 -n 'c (H$$_2$$C=CH$$_2$$)' -o ../FigureS6a.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/C2H4_ads.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a C3H6 -n 'd (H$$_2$$C=CH$$_2$$-CH$$_3$$)' -o ../FigureS6b.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/C3H6_ads.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a CO2 -n 'e (CO$$_2$$)' -o ../FigureS7a.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/CO2_ads.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a CH2O -n 'f (H$$_2$$C=O)' -o ../FigureS7b.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/CH2O_ads.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a C2H4O -n 'g (CH$$_3$$-CH=O)' -o ../FigureS8a.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/C2H4O_ads.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a THF -n 'h (THF)' -o ../FigureS8b.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/THF_ads.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a CH4O -n 'i (H$$_3$$C-OH)' -o ../FigureS9a.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/CH4O_ads.txt | tr '\n'  ',')"
	python ./plot_XPS_slabs_adsorbate_2.py -a THFn -n 'j (THF$$_n$$)' -o ../FigureS9b.pdf -x "$$(cat annotations/slabs_annotate.txt | tr '\n'  ','),$$(cat annotations/THFn_ads.txt | tr '\n'  ',')"

clean:
	rm ../FigureS*.pdf ../Figure{3,4,5,6,7a,7b}.pdf
