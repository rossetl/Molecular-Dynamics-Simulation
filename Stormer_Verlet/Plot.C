void Plot (double density, double energy)
{
	TCanvas* c1 = new TCanvas ("c1", "Pair Correlation Function", 800., 600.);
	TGraph* g = new TGraph ("ResultPCF.dat");
	g->SetLineColor(kRed);
	g->SetLineWidth(3);
	g->SetMarkerStyle(4);
	g->SetMarkerSize(.8);
	g->SetMarkerColor(kBlue);
	g->GetXaxis()->SetTitle("r");
	g->GetYaxis()->SetTitle("g(r)");
	string name = "Densita': " + to_string(density) + "       Energia: " + to_string(energy);
	g->SetTitle(name.c_str());
	g->Draw("ACE");
	c1->Draw();
}
