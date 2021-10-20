void PlotPCF (float density, float temperature)
{
	TCanvas* c1 = new TCanvas ("c1", "Pair Correlation Function", 800., 600.);
	TGraph* g = new TGraph ("ResultPCF.dat");
	g->SetLineColor(kRed);
	g->SetLineWidth(3);
	g->SetMarkerStyle(4);
	g->SetMarkerSize(.8);
	g->SetMarkerColor(kBlue);
	g->GetXaxis()->SetTitle("r/#sigma");
	g->GetYaxis()->SetTitle("g(r)");
	string name = "#rho: " + to_string(density) + "     T: " + to_string(temperature);
	g->SetTitle(name.c_str());
	g->Draw("ACE");
	c1->Draw();
}

void PlotEnergy()
{
	TCanvas* c1 = new TCanvas ("c1", "Energia", 800., 600.);
	TGraph* g = new TGraph ("Result.dat");
	g->SetLineColor(kBlue);
	g->SetLineWidth(3);
	g->SetMarkerStyle(4);
	g->SetMarkerSize(.8);
	g->GetXaxis()->SetTitle("step");
	g->GetYaxis()->SetTitle("E");
	g->SetTitle("Energia");
	g->Draw("ACE");
	c1->Draw();
}
