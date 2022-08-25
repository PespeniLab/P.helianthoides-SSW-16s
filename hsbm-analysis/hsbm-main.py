from seastar_hsbm_inference import seastar_hsbm

comm = seastar_hsbm()
comm.make_graph("level-7.csv")
comm.fit(deg_corr=False)
comm.save_model("hsbm-fit-non-dc")
comm.save_model_dict("fit-dict-non-dc")
comm.plot("plots/hsbm-fit-non-dc.svg")

comm = seastar_hsbm()
comm.make_graph("level-7.csv")
comm.fit(deg_corr=True)
comm.save_model("hsbm-fit")
comm.save_model_dict("fit-dict")
comm.plot("plots/hsbm-fit.svg")

