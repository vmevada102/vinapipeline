First time run following command to start conda

conda activate vina_pipeline


Start docking (First time)
./vina_pipeline.sh 

Resume
./vina_pipeline.sh --resume


After complition of the docking
python all_complexes_advanced_analysis.py 

Generate ranking
python analysis_pipeline.py

Final Screening
python final_parallel_screener_stable.py 
