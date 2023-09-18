for ont in mf bp cc; do
    echo $ont
    scp data/$ont/terms.pkl glogin.ibex.kaust.edu.sa:/ibex/user/kulmanm/deepgo-cafa5/data/$ont/
    scp data/$ont/interpros.pkl glogin.ibex.kaust.edu.sa:/ibex/user/kulmanm/deepgo-cafa5/data/$ont/
    scp data/$ont/*_data.pkl glogin.ibex.kaust.edu.sa:/ibex/user/kulmanm/deepgo-cafa5/data/$ont/
    scp data/$ont/*_data_diam.pkl glogin.ibex.kaust.edu.sa:/ibex/user/kulmanm/deepgo-cafa5/data/$ont/
done

