#include "cuts.h"

#include "TChain.h"
#include "TTreePlayer.h"

void makeCands() {
    //
    // Lepton+X
    //

    TChain *emu_chain = new TChain("tree");
    // emu_chain->Add("/store/disk00/jribnik/hunt/emu_baby/*.root");
    //emu_chain->Add("/store/disk00/jribnik/hunt/prompt_emu_baby/*.root");
    //emu_chain->Add("/nfs-3/userdata/yanjuntu/hunt/emu_baby/*.root");
    emu_chain->Add("emu_baby/*.root");
    TTreePlayer *emu_tp = (TTreePlayer*)emu_chain->GetPlayer();
    emu_tp->SetScanRedirect(kTRUE);
    //emu_tp->SetScanFileName("cands/wplusjetsv1_emu.raw.cands");
    //emu_chain->Scan("*", wplusjetsv1_emu, "colsize=50");
    emu_tp->SetScanFileName("cands/wprime100_emu.raw.cands");
    emu_chain->Scan("*", wprime100_emu, "colsize=50");
    //emu_tp->SetScanFileName("cands/semileptonictop_emu.raw.cands");
    //emu_chain->Scan("*", semileptonictop_emu, "colsize=50");
    //emu_tp->SetScanFileName("cands/semileptonictopv2_emu.raw.cands");
    //emu_chain->Scan("*", semileptonictopv2_emu, "colsize=50");
    emu_tp->SetScanFileName("cands/semileptonictopv3_emu.raw.cands");
    emu_chain->Scan("*", semileptonictopv3_emu, "colsize=50");
    //emu_tp->SetScanFileName("cands/highptmuon_emu.raw.cands");
    //emu_chain->Scan("*", highptmuon_emu, "colsize=50");
    emu_tp->SetScanFileName("cands/highptmuonv2_emu.raw.cands");
    emu_chain->Scan("*", highptmuonv2_emu, "colsize=50");
    emu_tp->SetScanFileName("cands/multimuon_emu.raw.cands");
    emu_chain->Scan("*", multimuon_emu, "colsize=50");

    //
    // Dilepton+X
    //

    TChain *dilep_chain = new TChain("tree");
    //dilep_chain->Add("/store/disk00/jribnik/hunt/dilep_baby/*.root");
    //dilep_chain->Add("/store/disk00/jribnik/hunt/prompt_dilep_baby/*.root");
    //dilep_chain->Add("/nfs-3/userdata/yanjuntu/hunt_old/dilep_baby/*.root");
    dilep_chain->Add("dilep_baby/*.root");
    TTreePlayer *dilep_tp = (TTreePlayer*)dilep_chain->GetPlayer();
    dilep_tp->SetScanRedirect(kTRUE);
    //dilep_tp->SetScanFileName("cands/inclusivez_dilep.raw.cands");
    //dilep_chain->Scan("*", inclusivez_dilep, "colsize=50");
    //dilep_tp->SetScanFileName("cands/zplusjetsv1_dilep.raw.cands");
    //dilep_chain->Scan("*", zplusjetsv1_dilep, "colsize=50");
    //dilep_tp->SetScanFileName("cands/zplusjetsv2_dilep.raw.cands");
    //dilep_chain->Scan("*", zplusjetsv2_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/zplusjetsv3_dilep.raw.cands");
    dilep_chain->Scan("*", zplusjetsv3_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/zplusbjetsv1_dilep.raw.cands");
    dilep_chain->Scan("*", zplusbjetsv1_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/zplusbjetsv2_dilep.raw.cands");
    dilep_chain->Scan("*", zplusbjetsv2_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/zmet30_dilep.raw.cands");
    dilep_chain->Scan("*", zmet30_dilep, "colsize=50");
    //dilep_tp->SetScanFileName("cands/zptz50_dilep.raw.cands");
    //dilep_chain->Scan("*", zptz50_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/zptz100_dilep.raw.cands");
    dilep_chain->Scan("*", zptz100_dilep, "colsize=50");
    //dilep_tp->SetScanFileName("cands/inclusivenonz_dilep.raw.cands");
    //dilep_chain->Scan("*", inclusivenonz_dilep, "colsize=50");
    //dilep_tp->SetScanFileName("cands/inclusivenonzhighmass_dilep.raw.cands");
    //dilep_chain->Scan("*", inclusivenonzhighmass_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/inclusivenonzhighmassv2_dilep.raw.cands");
    dilep_chain->Scan("*", inclusivenonzhighmassv2_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/inclusivenonzof_dilep.raw.cands");
    dilep_chain->Scan("*", inclusivenonzof_dilep, "colsize=50");
    //dilep_tp->SetScanFileName("cands/inclusivenonzjets_dilep.raw.cands");
    //dilep_chain->Scan("*", inclusivenonzjets_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/inclusivenonzjetsv2_dilep.raw.cands");
    dilep_chain->Scan("*", inclusivenonzjetsv2_dilep, "colsize=50");
    //dilep_tp->SetScanFileName("cands/inclusivenonzmet20_dilep.raw.cands");
    //dilep_chain->Scan("*", inclusivenonzmet20_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/inclusivenonzmet20nojets_dilep.raw.cands");
    dilep_chain->Scan("*", inclusivenonzmet20nojets_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/inclusivenonzmet40jets_dilep.raw.cands");
    dilep_chain->Scan("*", inclusivenonzmet40jets_dilep, "colsize=50");
    //dilep_tp->SetScanFileName("cands/inclusivenonzptz50_dilep.raw.cands");
    //dilep_chain->Scan("*", inclusivenonzptz50_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/inclusivenonzptz100_dilep.raw.cands");
    dilep_chain->Scan("*", inclusivenonzptz100_dilep, "colsize=50");
    //dilep_tp->SetScanFileName("cands/dileptonictop_dilep.raw.cands");
    //dilep_chain->Scan("*", dileptonictop_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/dileptonictopv2_dilep.raw.cands");
    dilep_chain->Scan("*", dileptonictopv2_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/dileptonictopv3_dilep.raw.cands");
    dilep_chain->Scan("*", dileptonictopv3_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/dileptonictopv4_dilep.raw.cands");
    dilep_chain->Scan("*", dileptonictopv4_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/samesigninclusive_dilep.raw.cands");
    dilep_chain->Scan("*", samesigninclusive_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/zplusworzprime_dilep.raw.cands");
    dilep_chain->Scan("*", zplusworzprime_dilep, "colsize=50");
    dilep_tp->SetScanFileName("cands/meff_dilep.raw.cands");
    dilep_chain->Scan("*", meff_dilep, "colsize=50");

    //
    // Trilepton+X
    //

    TChain *trilep_chain = new TChain("tree");
    //trilep_chain->Add("/store/disk00/jribnik/hunt/trilep_baby/*.root");
    //trilep_chain->Add("/store/disk00/jribnik/hunt/prompt_trilep_baby/*.root");
    //trilep_chain->Add("/nfs-3/userdata/yanjuntu/hunt_old/trilep_baby/*.root");
    trilep_chain->Add("trilep_baby/*.root");
    TTreePlayer *trilep_tp = (TTreePlayer*)trilep_chain->GetPlayer();
    trilep_tp->SetScanRedirect(kTRUE);
    trilep_tp->SetScanFileName("cands/inclusive_trilep.raw.cands");
    trilep_chain->Scan("*", inclusive_trilep, "colsize=50");
}
