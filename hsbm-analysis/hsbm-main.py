from seastar_hsbm_inference import seastar_hsbm

comm = seastar_hsbm()
comm.make_graph("level-7.csv")
comm.fit(deg_corr=False)
comm.save_model("hsbm-fit-non-dc2")
comm.save_model_dict("fit-dict-non-dc2")
comm.plot("hsbm-fit-non-dc2.pdf")

comm = seastar_hsbm()
comm.make_graph("level-7.csv")
comm.fit(deg_corr=True)
comm.save_model("hsbm-fit2")
comm.save_model_dict("fit-dict2")
comm.plot("hsbm-fit2.pdf")

# comm.load_model("hsbm-tmp-fit")

