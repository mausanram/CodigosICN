//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B02DetectorConstruction.hh"
#include "B02DetectorMessenger.hh"
#include "B02BarSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Paraboloid.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include "G4GenericMessenger.hh"
#include "G4VisAttributes.hh"
#include "G4VSensitiveDetector.hh"
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



B02DetectorConstruction::B02DetectorConstruction()

{
  DefineCommands();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02DetectorConstruction::~B02DetectorConstruction()
{
delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double B02DetectorConstruction::GetWorldFullLength() const
{return fWorldLength;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//added 


G4VPhysicalVolume* B02DetectorConstruction::Construct()
{ 
  //------------------------------------------------------ volumes
  //--------- Sizes of the principal geometrical components (solids)
  
  // new introductions to better comprehension
  DefineMaterials();
  
  // ============ Materials and Elements ======================== //
  // Air //
  auto Air = G4Material::GetMaterial("G4_AIR");
  auto LAr = new G4Material("LAr",   18.,  39.95*g/mole,  1.393*g/cm3);
  //LAr->GetIonisation()->SetBirksConstant(0.069*cm/MeV); //LarSoft documentation C. Zhang
  
  // Si //
  auto nistManager = G4NistManager::Instance();
  auto Si = nistManager->FindOrBuildMaterial("G4_Si");
  Si->GetIonisation()->SetBirksConstant(0.09*cm/MeV); // El valor de 0.1 parece ajustar bastante bien el espectro pero se necesita mas estudio
  
  //auto Si = new G4Material("Si",   14.,  28.0849*g/mole,  2.3396*g/cm3);
  //auto Si = G4Material::GetMaterial("G4_Si");
  //G4Element *elSi = new G4Element("Silicon", "Si", 14., 28.0849*g/mole);
  
  // Vacuum //
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  auto Vacuum = new G4Material("interGalactic", atomicNumber, massOfMole, density, kStateGas, temperature, pressure);
  
  // Steel (density ~8.0 g/cm³) //
  auto Steel = new G4Material("StainlessSteel", 8.0*g/cm3, 4);
  // Add elements with their fractional mass
  Steel->AddElement(elFe, 0.70);  // 70% Iron
  Steel->AddElement(elCr, 0.18);  // 18% Chromium
  Steel->AddElement(elNi, 0.10);  // 10% Nickel
  Steel->AddElement(elC,  0.02);  // 2% Carbon

  // Copper // 
  auto Cu = nistManager->FindOrBuildMaterial("G4_Cu");

  // Aluminum // 
  auto Al = nistManager->FindOrBuildMaterial("G4_Al");

  // AlN //
  G4Element* elAl = nistManager->FindOrBuildElement("Al");
  G4Element* elN = nistManager->FindOrBuildElement("N");
  G4Material* AlN = new G4Material("AluminumNitride", 3.26 * g/cm3, 2);
  AlN->AddElement(elAl, 1);
  AlN->AddElement(elN, 1);

  // ====================== // end materials
    
  fWorldLength = 150*cm;
  fRockLength = 0.80*fWorldLength;
  fRock2Length = 0.50*fWorldLength;
  fBarDiameter = 0.3*fRock2Length; 
     
  // World
  //------------------------------ 
  G4double HalfWorldLength = 1.*fWorldLength;
  //G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  //G4cout << "Computed tolerance = "
    //     << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
    //     << " mm" << G4endl;
  solidWorld= new G4Box("World",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  //logicWorld= new G4LogicalVolume( solidWorld, Vacuum, "World", 0, 0, 0);
  logicWorld= new G4LogicalVolume(solidWorld, Air, "World");
  physiWorld = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.), logicWorld, "World", 0, false, 0, true);   
  // =========================================== //
  
  // ============= Colors ================ //
  G4VisAttributes blue(G4Colour::Blue());
  G4VisAttributes cgray(G4Colour::Gray());
  G4VisAttributes green(G4Colour::Green());
  G4VisAttributes red(G4Colour::Red());
  G4VisAttributes yellow(G4Colour::Yellow());
  G4VisAttributes brown(G4Colour::Brown());
  G4VisAttributes white(G4Colour::White());
  
  G4double fact_transparent = 0.1;
  G4VisAttributes* semiTransparentYellow = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, fact_transparent));
  G4VisAttributes* semiTransparentLighBlue = new G4VisAttributes(G4Colour(173., 216., 230., fact_transparent));
  // ====================================== //

  G4double half_inch = 2.54/2; // cm

  // ================= SHIELDING ====================== //

  // ===== STEEL BOX ==== //
  G4bool flag_steelshielding = true;
  G4bool flag_coppershielding = true;

  G4double halfsize_box = 12.7; // cm
  G4double stbox_thickness = halfsize_box - 2.54/2; // cm

  auto box_st = new G4Box("box_st", halfsize_box*cm, halfsize_box*cm, halfsize_box*cm);
  // auto box_stLV = new G4LogicalVolume(box_st, Steel, "box_stLV");
  // //new G4PVPlacement(0, G4ThreeVector(0*cm, boxside/2-wth/2, 0*cm), yfLV1, "yf1", logicWorld, false, 0, fCheckOverlaps);

  G4double RadioCyl = halfsize_box- half_inch; 
  G4Tubs* cylinder1cut = new G4Tubs("cyl1", 0*cm, RadioCyl*cm, (halfsize_box+0.1)*cm,0,2*pi); // This function use ("name", in_radio, out_radio, high, angle_rad, angle_rad)
  G4SubtractionSolid* truebox_st = new G4SubtractionSolid("truebox_st", box_st, cylinder1cut, 0, G4ThreeVector((0)*cm, 0.*cm, (0)*cm));
  G4LogicalVolume* box_stLV = new G4LogicalVolume(truebox_st, Steel, "box_stLV");
  new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), box_stLV, "box_st", logicWorld, false, 0, true);

  auto vacuum_st = new G4Box("vacuum_st", stbox_thickness*cm, stbox_thickness*cm, stbox_thickness*cm);
  auto vacuum_stLV = new G4LogicalVolume(vacuum_st, Vacuum, "vacuum_stLV");
  new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), vacuum_stLV, "vacuum_stLV", box_stLV, false, 0, true);
  box_stLV->SetVisAttributes(blue);
  vacuum_stLV->SetVisAttributes(semiTransparentYellow);

  // Cylinders //
  // G4double RadioCyl = halfsize_box- half_inch; 
  auto cylinder1 = new G4Tubs("cyl1", 0*cm, RadioCyl*cm, (halfsize_box)*cm,0,2*pi); // This function use ("name", in_radio, out_radio, high, angle_rad, angle_rad)
  auto cyl1LV = new G4LogicalVolume(cylinder1, Vacuum, "cyl1LV");
  new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), cyl1LV, "cyl1LV", box_stLV, false, 0, true);
  cyl1LV->SetVisAttributes(semiTransparentYellow);

  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateX(90.*deg); 
  auto cylinder2 = new G4Tubs("cyl2", 0*cm, RadioCyl*cm, halfsize_box*cm,0,2*pi);
  auto cyl2LV = new G4LogicalVolume(cylinder2, Vacuum, "cyl2LV");
  new G4PVPlacement(rm, G4ThreeVector(0*cm, 0*cm, 0*cm), cyl2LV, "cyl2LV", box_stLV, false, 0, true);
  cyl2LV->SetVisAttributes(semiTransparentYellow);

  rm = new G4RotationMatrix();
  rm->rotateY(90.*deg);
  auto cylinder3 = new G4Tubs("cyl3", 0*cm, RadioCyl*cm, halfsize_box*cm,0,2*pi);
  auto cyl3LV = new G4LogicalVolume(cylinder3, Vacuum, "cyl3LV");
  new G4PVPlacement(rm, G4ThreeVector(0*cm, 0*cm, 0*cm), cyl3LV, "cyl3LV", box_stLV, false, 0, true);
  cyl3LV->SetVisAttributes(semiTransparentYellow);
  // ========= // end cylinders

  // Flaps //
  RadioCyl = halfsize_box;
  G4double flap_thickness = half_inch/2;

  rm = new G4RotationMatrix();
  rm->rotateY(90.*deg);
  auto flapx1 = new G4Tubs("flapx1", 0*cm, RadioCyl*cm, flap_thickness*cm,0,2*pi); // This function use ("name", in_radio, out_radio, high, angle_rad, angle_rad)
  auto flapx1LV = new G4LogicalVolume(flapx1, Steel, "flapx1LV");
  new G4PVPlacement(rm, G4ThreeVector((halfsize_box+flap_thickness)*cm, 0*cm, 0*cm), flapx1LV, "flapx1LV", box_stLV, false, 0, true);
  flapx1LV->SetVisAttributes(blue);

  auto flapx2 = new G4Tubs("flapx2", 0*cm, RadioCyl*cm, flap_thickness*cm,0,2*pi);
  auto flapx2LV = new G4LogicalVolume(flapx2, Steel, "flapx2LV");
  new G4PVPlacement(rm, G4ThreeVector(-(halfsize_box+flap_thickness)*cm, 0*cm, 0*cm), flapx2LV, "flapx2LV", box_stLV, false, 0, true);
  flapx2LV->SetVisAttributes(blue);

  //////////////// Comment all this block if you implement the cryocooler ///////////////
  // auto flapz1 = new G4Tubs("flapz1", 0*cm, RadioCyl*cm, flap_thickness*cm,0,2*pi); // This function use ("name", in_radio, out_radio, high, angle_rad, angle_rad)
  // auto flapz1LV = new G4LogicalVolume(flapz1, Steel, "flapz1LV");
  // new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, (halfsize_box+flap_thickness)*cm), flapz1LV, "flapz1LV", box_stLV, false, 0, true);
  // flapz1LV->SetVisAttributes(blue); 
  ///////////////////////////////////////////////////////////////////////////////////////

  auto flapz2 = new G4Tubs("flapz2", 0*cm, RadioCyl*cm, half_inch*cm,0,2*pi);
  auto flapz2LV = new G4LogicalVolume(flapz2, Steel, "flapz2LV");
  new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, -(halfsize_box+flap_thickness)*cm), flapz2LV, "flapz2LV", box_stLV, false, 0, true);
  flapz2LV->SetVisAttributes(blue);

  rm = new G4RotationMatrix();
  rm->rotateX(90.*deg);
  auto flapy1 = new G4Tubs("flapy1", 0*cm, RadioCyl*cm, flap_thickness*cm,0,2*pi); // This function use ("name", in_radio, out_radio, high, angle_rad, angle_rad)
  auto flapy1LV = new G4LogicalVolume(flapy1, Steel, "flapy1LV");
  new G4PVPlacement(rm, G4ThreeVector(0*cm, (halfsize_box+flap_thickness)*cm, 0*cm), flapy1LV, "flapy1LV", box_stLV, false, 0, true);
  flapy1LV->SetVisAttributes(blue);

  auto flapy2 = new G4Tubs("flapy2", 0*cm, RadioCyl*cm, flap_thickness*cm,0,2*pi);
  auto flapy2LV = new G4LogicalVolume(flapy2, Steel, "flapy2LV");
  new G4PVPlacement(rm, G4ThreeVector(0*cm, -(halfsize_box+flap_thickness)*cm, 0*cm), flapy2LV, "flapy2LV", box_stLV, false, 0, true);
  flapy2LV->SetVisAttributes(blue);
  // ====== END STEEL BOX === //


  // ====== COLD HEAD ========== //
  // Cylinder //
  RadioCyl = halfsize_box; 
  G4double half_high_cryoshield = 17.6/2 + half_inch/2; // cm

  G4Tubs* cryo_shield = new G4Tubs("cryo_shield", (RadioCyl - 2.54/2)*cm, RadioCyl*cm, (half_high_cryoshield)*cm,0,2*pi); // This function use ("name", in_radio, out_radio, high, angle_rad, angle_rad)
  G4Tubs* cryo_shieldcut = new G4Tubs("cryo_shieldcut", (RadioCyl - (half_inch)/2)*cm, (RadioCyl+0.1)*cm, (half_high_cryoshield-half_inch)*cm,0,2*pi);
  G4SubtractionSolid* truecryo_shield = new G4SubtractionSolid("truecryo_shield", cryo_shield, cryo_shieldcut, 0, G4ThreeVector((0)*cm, 0.*cm, (0)*cm));

  G4LogicalVolume* cryo_shieldLV = new G4LogicalVolume(truecryo_shield, Steel, "cryo_shieldLV");
  new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, (halfsize_box+half_high_cryoshield)*cm), cryo_shieldLV, "cryo_shieldLV", box_stLV, false, 0, true);
  cryo_shieldLV->SetVisAttributes(blue);

  G4Tubs* cryo_shieldfill = new G4Tubs("cryo_shieldfill", 0,(RadioCyl - half_inch-0.01)*cm, (half_high_cryoshield)*cm,0,2*pi);
  G4LogicalVolume* cryo_shieldfillLV = new G4LogicalVolume(cryo_shieldfill, Vacuum, "cryo_shieldfill");
  new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, (halfsize_box+half_high_cryoshield)*cm), cryo_shieldfillLV, "cryo_shieldfillLV", box_stLV, false, 0, true);
  cryo_shieldfillLV->SetVisAttributes(semiTransparentYellow);
  // End Cylinder //

  // Flange //
  G4double inRadio_flag = 17.4/2; // cm
  G4Tubs* flange = new G4Tubs("flange", (RadioCyl - inRadio_flag)*cm, RadioCyl*cm, (half_inch/2)*cm,0,2*pi); // This function use ("name", in_radio, out_radio, high, angle_rad, angle_rad)
  G4LogicalVolume* flangeLV = new G4LogicalVolume(flange, Steel, "flangeLV");
  new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, (halfsize_box+2*half_high_cryoshield+half_inch/2)*cm), flangeLV, "flangeLV", logicWorld, false, 0, true);
  flangeLV->SetVisAttributes(blue);

  G4Tubs* flangefill = new G4Tubs("flangefill", 0, (RadioCyl - inRadio_flag-0.01)*cm, (half_inch/2)*cm,0,2*pi);
  G4LogicalVolume* flangefillLV = new G4LogicalVolume(flangefill, Vacuum, "flangefill");
  new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, (halfsize_box+2*half_high_cryoshield+half_inch/2)*cm), flangefillLV, "flangefillLV", box_stLV, false, 0, true);
  flangefillLV->SetVisAttributes(semiTransparentYellow);
  // End flange

  // // == Coldhead == //
  // // Interface //
  // G4double Radio_interface = 12.34/2; // cm
  // G4double halfinterface_size = 5.39/2; // cm
  // G4Tubs* interface = new G4Tubs("interface", (Radio_interface - halfinterface_size)*cm, Radio_interface*cm, (half_inch/2)*cm,0,2*pi);
  // G4LogicalVolume* interfaceLV = new G4LogicalVolume(interface, Steel, "interfaceLV");
  // new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, (halfsize_box+2*half_high_cryoshield+1*half_inch+(half_inch/2))*cm), interfaceLV, "interfaceLV", logicWorld, false, 0, true);
  // interfaceLV->SetVisAttributes(blue);

  // G4Tubs* interface_fill = new G4Tubs("interface", 0, (Radio_interface - halfinterface_size)*cm, (half_inch/2)*cm,0,2*pi);
  // G4LogicalVolume* interface_fillLV = new G4LogicalVolume(interface_fill, Vacuum, "interface_fillLV");
  // new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, (halfsize_box+2*half_high_cryoshield+1*half_inch+(half_inch/2))*cm), interface_fillLV, "interface_fillLV", logicWorld, false, 0, true);
  // interface_fillLV->SetVisAttributes(semiTransparentYellow);
  // // End Interface //

  // // Head //
  // G4double halfhigh_head = 20.71/2; // cm
  // G4double halfsize_head = 8.8/2; // cm
  // G4Box* head = new G4Box("head", halfsize_head*cm, halfsize_head*cm, halfhigh_head*cm);

  // G4double thick_head = halfsize_head - 0.5; // cm
  // G4Box* head_cut = new G4Box("head_cut", (thick_head)*cm, (thick_head)*cm, (halfhigh_head+2)*cm);
  // G4SubtractionSolid* true_head = new G4SubtractionSolid("true_head", head, head_cut, 0, G4ThreeVector(0*cm, 0.*cm, (-thick_head-0.1)*cm));
  // G4LogicalVolume* headLV = new G4LogicalVolume(true_head, Al, "headLV");
  // new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, (halfsize_box+2*half_high_cryoshield+half_inch+(half_inch)+halfhigh_head)*cm), headLV, "headLV", logicWorld, false, 0, true);
  
  // G4Box* headfill = new G4Box("headfill", (thick_head-0.01)*cm, (thick_head-0.01)*cm, (halfhigh_head-2+0.7)*cm);
  // G4LogicalVolume* headfillLV = new G4LogicalVolume(headfill, Air, "headfillLV");
  // new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, (halfsize_box+2*half_high_cryoshield+5*half_inch+halfhigh_head/2)*cm), headfillLV, "headfillLV", logicWorld, false, 0, true);
  // headfillLV->SetVisAttributes(blue);
  // // End Head //

  // // // Tubs //
  // G4double Radio_tubs = 3.03/2; // cm
  // G4double thick_tubs = 0.5/2;
  // G4double halfhigh_tubs = 15.65/2; // cm
  // G4Tubs* tube1 = new G4Tubs("tube1", (Radio_tubs - thick_tubs)*cm, Radio_tubs*cm, (halfhigh_tubs+half_inch*4)*cm,0,2*pi);
  // G4LogicalVolume* tube1LV = new G4LogicalVolume(tube1, Al, "tube1LV");
  // new G4PVPlacement(0, G4ThreeVector((Radio_tubs)*cm, 0*cm, (halfsize_box+half_high_cryoshield+(half_inch))*cm), tube1LV, "tube1LV", box_stLV, false, 0, true);
  // tube1LV->SetVisAttributes(cgray);

  // G4Tubs* tube2 = new G4Tubs("tube2", (Radio_tubs - thick_tubs)*cm, Radio_tubs*cm, (halfhigh_tubs+half_inch*4)*cm,0,2*pi);
  // G4LogicalVolume* tube2LV = new G4LogicalVolume(tube2, Al, "tube2LV");
  // new G4PVPlacement(0, G4ThreeVector((-Radio_tubs)*cm, 0*cm, (halfsize_box+half_high_cryoshield+(half_inch))*cm), tube2LV, "tube2LV", box_stLV, false, 0, true);
  // tube2LV->SetVisAttributes(cgray);
  // // // End tubs //

  // // // Tubs Coldhead //
  // Radio_tubs = 3.208/2; // cm
  // thick_tubs = 0.5/2;
  // halfhigh_tubs = 5./2; // cm
  // G4Tubs* coldhead1 = new G4Tubs("coldhead1", (Radio_tubs - thick_tubs)*cm, Radio_tubs*cm, (halfhigh_tubs)*cm,0,2*pi);
  // G4LogicalVolume* coldhead1LV = new G4LogicalVolume(coldhead1, Cu, "coldhead1LV");
  // new G4PVPlacement(0, G4ThreeVector((-Radio_tubs)*cm, 0*cm, (3*halfhigh_tubs+half_inch/4 + 0.2)*cm), coldhead1LV, "coldhead1LV", box_stLV, false, 0, true);
  // coldhead1LV->SetVisAttributes(brown);

  // G4Tubs* coldhead2 = new G4Tubs("coldhead2", (Radio_tubs - thick_tubs)*cm, Radio_tubs*cm, (halfhigh_tubs)*cm,0,2*pi);
  // G4LogicalVolume* coldhead2LV = new G4LogicalVolume(coldhead2, Cu, "coldhead2LV");
  // new G4PVPlacement(0, G4ThreeVector((Radio_tubs)*cm, 0*cm, (3*halfhigh_tubs+half_inch/4 + 0.2)*cm), coldhead2LV, "coldhead2LV", box_stLV, false, 0, true);
  // coldhead2LV->SetVisAttributes(brown);
  // // // End tubs Colhead //

  // // // Base Coldhead //
  // G4double Radio_basecold = 7.06/2; // cm
  // G4double halfhigh_basecold = 0.77/2; // cm
  // G4Tubs* base_cold = new G4Tubs("base_cold", (0)*cm, Radio_basecold*cm, (halfhigh_basecold)*cm,0,2*pi);
  // G4LogicalVolume* base_coldLV = new G4LogicalVolume(base_cold, Cu, "base_coldLV");
  // new G4PVPlacement(0, G4ThreeVector((0)*cm, 0*cm, (2*halfhigh_tubs+0.1)*cm), base_coldLV, "base_coldLV", box_stLV, false, 0, true);
  // base_coldLV->SetVisAttributes(brown);
  // // ======  END COLD HEAD ===== //


  // ====== COPPER BASE ===== //
  rm = new G4RotationMatrix();
  rm->rotateX(180.*deg);
  G4double rot_value = -1.;
  G4double traslation = 4.6;

  // Base 1 //
  G4double halfx_size = 11.1/2; // cm
  G4double halfy_size = 10.2/2; // cm
  G4double halfthick_base = 0.2/2; // cm
  auto base_cu = new G4Box("base_cu", halfx_size*cm, halfy_size*cm, halfthick_base*cm);

  G4double cut_halfx_size = 7.722/2; // cm
  G4double cut_halfy_size = 8.992/2; // cm
  G4double cut_halfthick_base = 0.089/2; // cm
  auto cutbase_cu = new G4Box("cutbase_cu", cut_halfx_size*cm, cut_halfy_size*cm, cut_halfthick_base*cm);
  G4SubtractionSolid* cut1_baseCu = new G4SubtractionSolid("cut1_baseCu", base_cu, cutbase_cu, 0, G4ThreeVector((-0.604)*cm, 0.*cm, (cut_halfthick_base*2)*cm));

  G4double cut1_halfx_size = 5.461/2; // cm
  G4double cut1_halfy_size = 5.461/2; // cm
  G4double cut1_halfthick_base = 0.111/2; // cm
  G4Box* cut2base_cu = new G4Box("cut2base_cu", cut1_halfx_size*cm, cut1_halfy_size*cm, cut1_halfthick_base*2*cm);
  G4SubtractionSolid* true_baseCu = new G4SubtractionSolid("true_baseCu", cut1_baseCu, cut2base_cu, 0, G4ThreeVector((-0.604)*cm, 0.*cm, (0)*cm));

  G4LogicalVolume* true_baseCuLV = new G4LogicalVolume(true_baseCu, Cu, "true_baseCu");
  new G4PVPlacement(rm, G4ThreeVector(0*cm, 0*cm, (traslation)*cm), true_baseCuLV, "true_baseCu", box_stLV, false, 0, true);
  true_baseCuLV->SetVisAttributes(brown);
  // End Base 1 //

  // Substrate //
  G4Box* substrate_AlN = new G4Box("substrate_AlN", cut_halfx_size*cm, cut_halfy_size*cm, (cut_halfthick_base*0.6)*cm);
  G4LogicalVolume* substrate_AlNLV = new G4LogicalVolume(substrate_AlN, AlN, "substrate_AlN");
  new G4PVPlacement(rm, G4ThreeVector((-0.604)*cm, 0.*cm, (traslation+ rot_value*(halfthick_base/2 + cut_halfthick_base/2))*cm), substrate_AlNLV, "substrate_AlN", box_stLV, false, 0, true);
  substrate_AlNLV->SetVisAttributes(cgray);
  // End Substrate //

  // Connector //
  G4double halfx_conn1 = 0.508/2; // cm
  G4double halfy_conn1 = 11.684/2; // cm
  G4double halfz_conn1 = 0.279/2; // cm

  G4Box*  Conn1 = new G4Box("Conn1", halfx_conn1*cm, halfy_conn1*cm, halfz_conn1*cm);
  G4LogicalVolume* Conn1LV = new G4LogicalVolume(Conn1, Cu, "Conn1LV");
  new G4PVPlacement(rm, G4ThreeVector((halfy_size+halfx_conn1*2.77+0.001)*cm, 0*cm, (traslation+rot_value*(-0.037))*cm), Conn1LV, "Conn1LV", box_stLV, false, 0, true);
  Conn1LV->SetVisAttributes(brown);

  G4double halfx_conn2 = 1.143/2; // cm
  G4double halfy_conn2 = 11.684/2; // cm
  G4double halfz_conn2 = 0.254/2; // cm
  G4Box*  Conn2 = new G4Box("Conn2", halfx_conn2*cm, halfy_conn2*cm, halfz_conn2*cm);

  G4double halfx_cut = 1.143/2; // cm
  G4double halfy_cut = 2.875/2; // cm
  G4double halfz_cut = 0.1524/2; // cm
  G4Box* CutConn2 = new G4Box("CutConn2", (halfx_cut+0.1)*cm, halfy_cut*cm, halfz_cut*cm);
  G4SubtractionSolid* true_Conn2 = new G4SubtractionSolid("true_Conn2", Conn2, CutConn2, 0, G4ThreeVector((0)*cm, 0.*cm, (halfz_cut*2)*cm));

  G4LogicalVolume* Conn2LV = new G4LogicalVolume(true_Conn2, Cu, "Conn2LV");
  new G4PVPlacement(rm, G4ThreeVector((halfy_size+halfx_conn1+0.136)*cm, 0*cm, (traslation+rot_value*(halfz_conn1+halfz_conn2-0.038))*cm), Conn2LV, "Conn2LV", box_stLV, false, 0, true);
  Conn2LV->SetVisAttributes(brown);
  // End Connector


  // Lid //
  G4double CCD_ZLength = 0.0725; // cm
  G4double halfx_lid = 8.573/2; // cm
  G4double halfy_lid = 10.16/2; // cm
  G4double halfz_lid = 0.167/2; // cm
  G4Box*  Lid = new G4Box("Lid", halfx_lid*cm, halfy_lid*cm, halfz_lid*cm);

  G4double halfx_lidcut = 6.9596/2; // cm
  G4double halfy_lidcut = 6.731/2; // cm
  G4double halfz_lidcut = halfz_lid; // cm
  G4Box* cutLid = new G4Box("cutLid", halfx_lidcut*cm, halfy_lidcut*cm, (halfz_lidcut+0.1)*cm);
  G4SubtractionSolid* true_Lid = new G4SubtractionSolid("true_Lid", Lid, cutLid, 0, G4ThreeVector((0)*cm, 0.*cm, (0)*cm));

  G4LogicalVolume* LidLV = new G4LogicalVolume(true_Lid, Cu, "LidLV");
  new G4PVPlacement(rm, G4ThreeVector((-0.5)*cm, 0*cm, (traslation+rot_value*(halfthick_base+halfz_lidcut))*cm), LidLV, "LidLV", box_stLV, false, 0, true);
  LidLV->SetVisAttributes(brown);

  // ===== END COPPER BASE ===== //


  // =============== Constructor of CCD (non active volume) ===================== //
  G4double pixel_size = 0.0015; // cm
  G4double XLength = 250 * pixel_size; // cm
  G4double YLength = 529 * pixel_size; // cm (Debe ser la dimensión con mayor tamaño)
  G4double ZLength = 0.0725; // cm
  // ============================================================================ //

  // ================== CCD ================================= //
  halfthick_base = 0.2/2;
  //Sibox = new G4Box("ccd", HalfWorldLength, HalfWorldLength, HalfWorldLength);
  Sibox = new G4Box("CCD", 0.5*XLength*cm, 0.5*YLength*cm, 0.5*ZLength*cm);
  SiLogic = new G4LogicalVolume(Sibox, Si, "CCD", 0, 0, 0);
  new G4PVPlacement(rm, G4ThreeVector(0., 0., (traslation+rot_value*(halfthick_base+ZLength/2))*cm), SiLogic, "CCD", box_stLV, false, 0,true);
  SiLogic->SetVisAttributes(white);
  fSiLogic = SiLogic;
  // ======================================================== //


  // =================== Sensitive detectors ================= //	
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String barSDname = "CCD/SD";
  B02BarSD* aBarSD = new B02BarSD( barSDname );
  SDman->AddNewDetector( aBarSD );
  SetSensitiveDetector("CCD", aBarSD, true);
  G4VisAttributes* BarVisAtt   = new G4VisAttributes (G4Colour(1,1,1));
  BarVisAtt->SetForceAuxEdgeVisible( true );
  SiLogic->SetVisAttributes(BarVisAtt);
  
 return physiWorld;
 
}


void B02DetectorConstruction::DefineMaterials()
{
 //------------------------------------------------------ materials
  
  
 //Some elements
    elFe = new G4Element("Iron", "Fe", 26., 55.85*g/mole);
    elCr = new G4Element("Chromium", "Cr", 24., 51.9961*g/mole);
    elNi = new G4Element("Nickel", "Ni", 28., 58.69*g/mole);
    elC  = new G4Element("Carbon", "C", 6., 12.01*g/mole);

//Some Nist materials
  
   auto nistManager = G4NistManager::Instance();

  // Air
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  //H2O
  nistManager->FindOrBuildMaterial("G4_WATER");
  
  //concrete 
  nistManager->FindOrBuildMaterial("G4_CONCRETE");
  
  //polyethylene
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  
  
  //G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
 
}

void B02DetectorConstruction::SetHigh(G4double val)
{
  if (!Sibox) {
      G4cerr << "Detector has not yet been constructed." << G4endl;
      return;
  }
  // zhigh = val;
  // SolidLAr-> SetZHalfLength(0.5*zhigh);
  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void B02DetectorConstruction::DefineCommands()
{
  
  // detectorMessenger = new B02DetectorMessenger(this);
  fMessenger = new G4GenericMessenger(this, "/detector/", "B02DetectorConstruction");
  
  fMessenger->DeclareMethodWithUnit("zhigh", "cm",  &B02DetectorConstruction::SetHigh, "high cylinder");
  
}


