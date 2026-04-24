{ 
    TFile f("combined_analysed.root");
    const char* names[] = {"h1_califa_gamma_corr_25F_p2p_22O_allmult_0","h1_califa_gamma_corr_25F_p2p_22O_mult_exact_0","h1_califa_gamma_corr_25F_p2p_22O_mult_leq_0","h1_califa_gamma_corr_25F_p2p_24O_0","h1_califa_gamma_corr_25F_p2p_24O_mult_exact_0","h1_califa_gamma_corr_25F_p2p_24O_mult_leq_0"};
    for(size_t i=0;i<sizeof(names)/sizeof(names[0]);++i){
        TObject* o = f.Get(names[i]);
        printf("%s %s\n", names[i], o?"FOUND":"MISSING");
    }
    f.Close();
}
