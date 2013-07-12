{
gROOT->ProcessLine(".L examples/AfbFinal2DUnfoldttrapidity2.cxx");
AfbUnfoldExample(1.,2.,1.,1.,1.,1.);
AfbUnfoldExample(1.,0.,1.,1.,1.,1.);
AfbUnfoldExample(1.,1.,2.,1.,1.,1.);
AfbUnfoldExample(1.,1.,0.,1.,1.,1.);
AfbUnfoldExample(1.,1.,1.,2.,1.,1.);
AfbUnfoldExample(1.,1.,1.,0.,1.,1.);
AfbUnfoldExample(1.,1.,1.,1.,1.5,1.);
AfbUnfoldExample(1.,1.,1.,1.,0.5,1.);
AfbUnfoldExample(1.,1.,1.,1.,1.,1.5);
AfbUnfoldExample(1.,1.,1.,1.,1.,0.5);
AfbUnfoldExample(1.,1.,1.,1.,1.,1.);
}
