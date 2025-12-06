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
#include "G4UnionSolid.hh"

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

  // STEEL BOX //
  G4double half_inch = 2.54/2 * cm;       // 1.27 cm
  G4double halfsize_box = 12.7 * cm;      
  G4double thick_wall = half_inch;
  G4double inner_halfsize = halfsize_box - thick_wall; 
  G4double RadioCyl = halfsize_box - half_inch; 

  G4Box* solid_Box_Outer = new G4Box("Box_Outer", halfsize_box, halfsize_box, halfsize_box);
  G4Box* solid_Box_Inner = new G4Box("Box_Inner", inner_halfsize, inner_halfsize, inner_halfsize);

  G4Tubs* solid_Drill = new G4Tubs("Drill", 0., RadioCyl, halfsize_box * 1.1, 0., 2*pi);

  // --- Rotations for the drill
  G4RotationMatrix* rotX = new G4RotationMatrix(); rotX->rotateX(90.*deg);
  G4RotationMatrix* rotY = new G4RotationMatrix(); rotY->rotateY(90.*deg);

  // --- Substractions 
  /// Box
  G4SubtractionSolid* steel_Step1 = new G4SubtractionSolid("Steel_Step1", solid_Box_Outer, solid_Box_Inner);

  /// Size Z
  G4SubtractionSolid* steel_Step2 = new G4SubtractionSolid("Steel_Step2", steel_Step1, solid_Drill);

  /// Size X
  G4SubtractionSolid* steel_Step3 = new G4SubtractionSolid("Steel_Step3", steel_Step2, solid_Drill, rotX, G4ThreeVector(0,0,0));

  /// Size Y
  G4SubtractionSolid* solid_FinalSteel = new G4SubtractionSolid("Steel_Final", steel_Step3, solid_Drill, rotY, G4ThreeVector(0,0,0));

  // --- End Substractions


  // We place the box IN THE WORLD
  G4LogicalVolume* logicSteelBox = new G4LogicalVolume(solid_FinalSteel, Steel, "Logic_SteelBox");
  new G4PVPlacement(0, G4ThreeVector(0,0,0), logicSteelBox, "Phys_SteelBox", logicWorld, false, 0, true);
  G4VisAttributes* attSteel = new G4VisAttributes(G4Colour::Blue());
  attSteel->SetForceAuxEdgeVisible(true);
  logicSteelBox->SetVisAttributes(attSteel);

  // --- Fill Construction
  /// Cylider
  G4Tubs* solid_VacCyl = new G4Tubs("VacCyl", 0., RadioCyl, halfsize_box, 0., 2*pi);

  /// Vacuum box + Zcyl
  G4UnionSolid* vac_Step1 = new G4UnionSolid("Vac_Step1", solid_Box_Inner, solid_VacCyl);

  /// + Ycyl
  G4UnionSolid* vac_Step2 = new G4UnionSolid("Vac_Step2", vac_Step1, solid_VacCyl, rotX, G4ThreeVector(0,0,0));

  /// + Xcyl
  G4UnionSolid* solid_FinalVacuum = new G4UnionSolid("Vacuum_Final", vac_Step2, solid_VacCyl, rotY, G4ThreeVector(0,0,0));

  // --- END Fill Construction

  G4LogicalVolume* logicVacuumInterior = new G4LogicalVolume(solid_FinalVacuum, Vacuum, "Logic_VacuumInterior");
  new G4PVPlacement(0, G4ThreeVector(0,0,0), logicVacuumInterior, "Phys_VacuumInterior", logicWorld, false, 0, true);
  G4VisAttributes* attVac = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.3));
  logicVacuumInterior->SetVisAttributes(attVac);

  // We save the pointer into the fill box, this will be the mother of other objects
  G4LogicalVolume* vacuum_stLV = logicVacuumInterior;

  // END STEEL BOX //

  // ========================================================
  //                       COLD HEAD
  // ========================================================
  G4double RadioCyl_Tower = halfsize_box; 
  G4double z_current_pos = halfsize_box;

  G4double Radio_tubs = 3.03/2 * cm;
  G4double thick_tubs = 0.5/2 * cm;
  G4double halfhigh_tubs = 15.65/2 * cm;
  
  G4Tubs* solid_TubeHole = new G4Tubs("TubeHole", 0., Radio_tubs + 0.001*cm, halfhigh_tubs * 2, 0., 2*pi);
  G4Tubs* solid_AlTube   = new G4Tubs("AlTube", Radio_tubs - thick_tubs, Radio_tubs, halfhigh_tubs, 0., 2*pi);


  // --------------------------------------------------------
  // A) CRYO SHIELD
  // --------------------------------------------------------
  G4double half_h_shield = 17.6/2 * cm + half_inch/2;
  z_current_pos += half_h_shield;


  G4Tubs* solid_Shield_St = new G4Tubs("Shield_St", RadioCyl_Tower - half_inch, RadioCyl_Tower, half_h_shield, 0, 2*pi);
  G4LogicalVolume* logicShieldSt = new G4LogicalVolume(solid_Shield_St, Steel, "Logic_Shield_St");
  new G4PVPlacement(0, G4ThreeVector(0,0, z_current_pos), logicShieldSt, "Phys_Shield_St", logicWorld, false, 0, true);

  // 2. Vacío Interior (CON AGUJEROS)
  // Cilindro de vacío bruto
  G4Tubs* solid_Shield_Vac_Raw = new G4Tubs("Shield_Vac_Raw", 0., RadioCyl_Tower - half_inch, half_h_shield, 0, 2*pi);
  

  // Restamos los agujeros para los tubos (ubicados en +/- Radio_tubs en X)
  // Importante: El agujero es largo, así que centrado en Z corta todo el cilindro
  G4SubtractionSolid* vac_Sh_Step1 = new G4SubtractionSolid("Vac_Sh_S1", solid_Shield_Vac_Raw, solid_TubeHole, 0, G4ThreeVector(Radio_tubs, 0, 0));
  G4SubtractionSolid* solid_Shield_Vac = new G4SubtractionSolid("Vac_Sh_Fin", vac_Sh_Step1, solid_TubeHole, 0, G4ThreeVector(-Radio_tubs, 0, 0));

  G4LogicalVolume* logicShieldVac = new G4LogicalVolume(solid_Shield_Vac, Vacuum, "Logic_Shield_Vac");
  new G4PVPlacement(0, G4ThreeVector(0,0, z_current_pos), logicShieldVac, "Phys_Shield_Vac", logicWorld, false, 0, true);

  logicShieldSt->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
  logicShieldVac->SetVisAttributes(G4VisAttributes(G4Colour(1., 1., 0., 0.3))); // Amarillo transp.

  z_current_pos += half_h_shield; // Movemos cursor al final del Shield


  // --------------------------------------------------------
  // B) FLANGE
  // --------------------------------------------------------
  G4double half_h_flange = half_inch / 2;
  G4double inRadio_flange = 17.4/2 * cm;
  z_current_pos += half_h_flange; 

  G4Tubs* solid_Flange_St = new G4Tubs("Flange_St", RadioCyl_Tower - inRadio_flange, RadioCyl_Tower, half_h_flange, 0, 2*pi);
  G4LogicalVolume* logicFlangeSt = new G4LogicalVolume(solid_Flange_St, Steel, "Logic_Flange_St");
  new G4PVPlacement(0, G4ThreeVector(0,0, z_current_pos), logicFlangeSt, "Phys_Flange_St", logicWorld, false, 0, true);

  G4Tubs* solid_Flange_Vac_Raw = new G4Tubs("Flange_Vac_Raw", 0., RadioCyl_Tower - inRadio_flange, half_h_flange, 0, 2*pi);  
  G4SubtractionSolid* vac_Fl_Step1 = new G4SubtractionSolid("Vac_Fl_S1", solid_Flange_Vac_Raw, solid_TubeHole, 0, G4ThreeVector(Radio_tubs, 0, 0));
  G4SubtractionSolid* solid_Flange_Vac = new G4SubtractionSolid("Vac_Fl_Fin", vac_Fl_Step1, solid_TubeHole, 0, G4ThreeVector(-Radio_tubs, 0, 0));

  G4LogicalVolume* logicFlangeVac = new G4LogicalVolume(solid_Flange_Vac, Vacuum, "Logic_Flange_Vac");
  new G4PVPlacement(0, G4ThreeVector(0,0, z_current_pos), logicFlangeVac, "Phys_Flange_Vac", logicWorld, false, 0, true);

  logicFlangeSt->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
  logicFlangeVac->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.3)));

  z_current_pos += half_h_flange;


  // --------------------------------------------------------
  // C) INTERFACE
  // --------------------------------------------------------
  G4double half_h_inter = half_inch / 2;
  G4double Radio_interface = 12.34/2 * cm;
  G4double halfinterface_size = 5.39/2 * cm;
  // z_current_pos += half_inch; 
  z_current_pos += half_h_inter;

  G4Tubs* solid_Inter_St = new G4Tubs("Inter_St", Radio_interface - halfinterface_size, Radio_interface, half_h_inter, 0, 2*pi);
  G4LogicalVolume* logicInterSt = new G4LogicalVolume(solid_Inter_St, Steel, "Logic_Inter_St");
  new G4PVPlacement(0, G4ThreeVector(0,0, z_current_pos), logicInterSt, "Phys_Inter_St", logicWorld, false, 0, true);

  G4Tubs* solid_Inter_Vac_Raw = new G4Tubs("Inter_Vac_Raw", 0., Radio_interface - halfinterface_size, half_h_inter, 0, 2*pi);  
  G4SubtractionSolid* vac_In_Step1 = new G4SubtractionSolid("Vac_In_S1", solid_Inter_Vac_Raw, solid_TubeHole, 0, G4ThreeVector(Radio_tubs, 0, 0));
  G4SubtractionSolid* solid_Inter_Vac = new G4SubtractionSolid("Vac_In_Fin", vac_In_Step1, solid_TubeHole, 0, G4ThreeVector(-Radio_tubs, 0, 0));

  G4LogicalVolume* logicInterVac = new G4LogicalVolume(solid_Inter_Vac, Vacuum, "Logic_Inter_Vac");
  new G4PVPlacement(0, G4ThreeVector(0,0, z_current_pos), logicInterVac, "Phys_Inter_Vac", logicWorld, false, 0, true);

  logicInterSt->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
  logicInterVac->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.3)));


  // --------------------------------------------------------
  // Aluminum tubs 
  // --------------------------------------------------------
  G4double z_tubes_pos = halfsize_box + half_h_shield + half_inch;

  G4LogicalVolume* logicAlTube = new G4LogicalVolume(solid_AlTube, Al, "Logic_AlTube");
  logicAlTube->SetVisAttributes(G4VisAttributes(G4Colour::Gray()));

  // Tubo 1 (+X)
  new G4PVPlacement(0, G4ThreeVector(Radio_tubs, 0, z_tubes_pos), logicAlTube, "Phys_Tube1", logicWorld, false, 0, true);
  // Tubo 2 (-X)
  new G4PVPlacement(0, G4ThreeVector(-Radio_tubs, 0, z_tubes_pos), logicAlTube, "Phys_Tube2", logicWorld, false, 1, true);


  // --------------------------------------------------------
  // Aluminum HEAD
  // --------------------------------------------------------
  G4double halfhigh_head = 20.71/2 * cm;
  G4double halfsize_head = 8.8/2 * cm;
  // Posicion Z: justo encima de la interfaz
  z_current_pos += half_h_inter + halfhigh_head; 

  // Caja exterior
  G4Box* solid_Head_Out = new G4Box("Head_Out", halfsize_head, halfsize_head, halfhigh_head);
  // Hueco interior
  G4double thick_head = halfsize_head - 0.5*cm; 
  G4Box* solid_Head_In  = new G4Box("Head_In", thick_head, thick_head, halfhigh_head); // Cortamos hasta abajo

  G4SubtractionSolid* solid_Head_Final = new G4SubtractionSolid("Head_Final", solid_Head_Out, solid_Head_In, 0, G4ThreeVector(0,0, -0.6*cm)); // Ajuste fino de fondo

  G4LogicalVolume* logicHead = new G4LogicalVolume(solid_Head_Final, Al, "Logic_Head");
  new G4PVPlacement(0, G4ThreeVector(0,0, z_current_pos), logicHead, "Phys_Head", logicWorld, false, 0, true);
  logicHead->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));
  // ==========================
  // END COLD HEAD
  // ==========================


  // ==========================
  // COLD HEAD INSIDE THE VACUUM
  // ==========================

  G4double Radio_cold_tubs = 3.208/2 * cm;
  G4double halfhigh_cold_tubs = 5.0/2 * cm;

  G4Tubs* solid_ColdExt = new G4Tubs("ColdExtension", 0., Radio_cold_tubs, halfhigh_cold_tubs, 0., 2*pi);
  G4LogicalVolume* logicColdExt = new G4LogicalVolume(solid_ColdExt, Cu, "Logic_ColdExt");
  logicColdExt->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));

  G4double z_pos_cold = inner_halfsize - halfhigh_cold_tubs - 0.05*cm; 
  // G4double posX_tubes = 3.03/2 * cm; 
  G4double posX_tubes = Radio_cold_tubs;

  new G4PVPlacement(0, G4ThreeVector(posX_tubes, 0, z_pos_cold), logicColdExt, "Phys_ColdExt1", vacuum_stLV, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(-posX_tubes, 0, z_pos_cold), logicColdExt, "Phys_ColdExt2", vacuum_stLV, false, 1, true);


  // CU BASE 
  G4double Radio_basecold = 7.06/2 * cm;
  G4double halfhigh_basecold = 0.77/2 * cm;

  G4Tubs* solid_BaseDisk = new G4Tubs("BaseDisk", 0., Radio_basecold, halfhigh_basecold, 0., 2*pi);
  G4LogicalVolume* logicBaseDisk = new G4LogicalVolume(solid_BaseDisk, Cu, "Logic_BaseDisk");
  logicBaseDisk->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));

  G4double z_pos_disk = z_pos_cold - halfhigh_cold_tubs - halfhigh_basecold;
  new G4PVPlacement(0, G4ThreeVector(0, 0, z_pos_disk), logicBaseDisk, "Phys_BaseDisk", vacuum_stLV, false, 0, true);

  // ==========================
  // END COLD HEAD INSIDE THE VACUUM
  // ==========================

  // ========================================================
  // CU BASE  + CCD
  // ========================================================
  G4double z_anchor_point = z_pos_disk - halfhigh_basecold; 
  G4double custom_Z_shift = -5.0 * cm; // if it's equal to 0 then the base and CCD will be close to the base cold head
  G4double z_start_base = z_anchor_point + custom_Z_shift;

  
  /// Cu Base 
  G4double base_x = 11.1/2 * cm;
  G4double base_y = 10.2/2 * cm;
  G4double base_z = 0.2/2 * cm; // Grosor

  G4RotationMatrix* rotFlip = new G4RotationMatrix();
  rotFlip->rotateX(180.*deg);

  // Base 1
  G4Box* solid_BaseMain = new G4Box("BaseMain", base_x, base_y, base_z);

  G4double cut_x = 7.722/2 * cm;
  G4double cut_y = 8.992/2 * cm;
  G4double cut_z = 0.089/2 * cm; // Profundidad corte
  G4Box* solid_CutPocket = new G4Box("CutPocket", cut_x, cut_y, cut_z * 2); // x2 para asegurar corte
  G4SubtractionSolid* solid_Base_Step1 = new G4SubtractionSolid("Base_S1", solid_BaseMain, solid_CutPocket, 0, G4ThreeVector(0, 0, base_z));

  G4double hole_x = 5.461/2 * cm;
  G4Box* solid_CenterHole = new G4Box("CenterHole", hole_x, hole_x, base_z * 2.1);
  G4SubtractionSolid* solid_Base_Final = new G4SubtractionSolid("Base_Final", solid_Base_Step1, solid_CenterHole, 0, G4ThreeVector(0,0,0));

  G4LogicalVolume* logicBaseCu = new G4LogicalVolume(solid_Base_Final, Cu, "Logic_BaseCu");
  logicBaseCu->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));

  // Posición de la Base:
  // Como rotamos 180 grados, el centro Z local se invierte.
  // Queremos que la cara "trasera" (que ahora está arriba) toque el disco.
  // Z_center = z_start_base - base_z (mitad grosor)
  G4double z_center_base = z_start_base - base_z;

  new G4PVPlacement(rotFlip, G4ThreeVector(0, 0, z_center_base), logicBaseCu, "Phys_BaseCu", vacuum_stLV, false, 0, true);


  // Substrate (AlN)
  G4Box* solid_Substrate = new G4Box("Substrate", cut_x, cut_y, cut_z); // Mismo tamaño que el corte
  G4LogicalVolume* logicSubstrate = new G4LogicalVolume(solid_Substrate, AlN, "Logic_Substrate");
  logicSubstrate->SetVisAttributes(G4VisAttributes(G4Colour::Gray()));

  G4double z_substrate_pos = z_center_base - (base_z - cut_z); 
  new G4PVPlacement(rotFlip, G4ThreeVector(0, 0, z_substrate_pos), logicSubstrate, "Phys_Substrate", vacuum_stLV, false, 0, true);


  // CCD
  G4double pixel_size = 0.0015 * cm; // 15 micras
  G4double ccd_x = 250 * pixel_size / 2;
  G4double ccd_y = 529 * pixel_size / 2;
  G4double ccd_z = 0.0725 / 2 * cm;

  G4Box* solid_CCD = new G4Box("CCD", ccd_x, ccd_y, ccd_z);
  G4LogicalVolume* logicCCD = new G4LogicalVolume(solid_CCD, Si, "Logic_CCD");
  logicCCD->SetVisAttributes(G4VisAttributes(G4Colour::White())); // Blanco brillante

  G4double z_ccd_pos = z_substrate_pos - cut_z - ccd_z;
  new G4PVPlacement(rotFlip, G4ThreeVector(0, 0, z_ccd_pos), logicCCD, "Phys_CCD", vacuum_stLV, false, 0, true);


  // --- F) Sensitive Detector (SD) ---  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String barSDname = "CCD/SD";
  B02BarSD* aBarSD = (B02BarSD*)SDman->FindSensitiveDetector(barSDname, false);
  if (!aBarSD) {
      aBarSD = new B02BarSD(barSDname);
      SDman->AddNewDetector(aBarSD);
  }
  logicCCD->SetSensitiveDetector(aBarSD);

  fSiLogic = logicCCD;


  // ========================================================
  // DETALLES INTERNOS (CONECTORES Y TAPA)
  // ========================================================

  // CONECTOR 1 ---
  G4double halfx_conn1 = 0.508/2 * cm;
  G4double halfy_conn1 = 11.684/2 * cm;
  G4double halfz_conn1 = 0.279/2 * cm;

  G4Box* solid_Conn1 = new G4Box("Conn1", halfx_conn1, halfy_conn1, halfz_conn1);
  G4LogicalVolume* logicConn1 = new G4LogicalVolume(solid_Conn1, Cu, "Logic_Conn1");
  logicConn1->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));

  G4double x_pos_conn1 = base_y + halfx_conn1*2.77 + 0.001*cm; 
  G4double z_pos_conn1 = z_center_base - (-0.037*cm); 
  new G4PVPlacement(rotFlip, G4ThreeVector(x_pos_conn1, 0, z_pos_conn1), logicConn1, "Phys_Conn1", vacuum_stLV, false, 0, true);


  // CONECTOR 2 ---
  G4double halfx_conn2 = 1.143/2 * cm;
  G4double halfy_conn2 = 11.684/2 * cm;
  G4double halfz_conn2 = 0.254/2 * cm;


  G4Box* solid_Conn2_Raw = new G4Box("Conn2_Raw", halfx_conn2, halfy_conn2, halfz_conn2);
  
  // Corte del conector
  G4double halfx_cut_c2 = 1.143/2 * cm;
  G4double halfy_cut_c2 = 2.875/2 * cm;
  G4double halfz_cut_c2 = 0.1524/2 * cm;
  G4Box* solid_CutConn2 = new G4Box("CutConn2", halfx_cut_c2+0.1*cm, halfy_cut_c2, halfz_cut_c2);

  G4SubtractionSolid* solid_Conn2 = new G4SubtractionSolid("Conn2_Final", solid_Conn2_Raw, solid_CutConn2, 0, G4ThreeVector(0, 0, halfz_cut_c2*2));

  G4LogicalVolume* logicConn2 = new G4LogicalVolume(solid_Conn2, Cu, "Logic_Conn2");
  logicConn2->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));

  G4double gap_safety = 0.002 * cm; // 20 micras de separación
  G4double x_pos_conn2 = base_y + halfx_conn1 + 0.136*cm; 
  G4double z_pos_conn2 = (z_pos_conn1 - halfz_conn1) - halfz_conn2 - gap_safety;

  new G4PVPlacement(rotFlip, G4ThreeVector(x_pos_conn2, 0, z_pos_conn2), logicConn2, "Phys_Conn2", vacuum_stLV, false, 0, true);


  // LID ---
  G4double halfx_lid = 8.573/2 * cm;
  G4double halfy_lid = 10.16/2 * cm;
  G4double halfz_lid = 0.167/2 * cm;

  G4Box* solid_Lid_Raw = new G4Box("Lid_Raw", halfx_lid, halfy_lid, halfz_lid);

  // Cut
  G4double halfx_lidcut = 6.9596/2 * cm;
  G4double halfy_lidcut = 6.731/2 * cm;
  G4Box* solid_LidCut = new G4Box("LidCut", halfx_lidcut, halfy_lidcut, halfz_lid+0.1*cm);
  G4SubtractionSolid* solid_Lid = new G4SubtractionSolid("Lid_Final", solid_Lid_Raw, solid_LidCut);
  
  G4LogicalVolume* logicLid = new G4LogicalVolume(solid_Lid, Cu, "Logic_Lid");
  logicLid->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));

  G4double z_pos_lid = z_center_base - (base_z + halfz_lid); 
  new G4PVPlacement(rotFlip, G4ThreeVector(-0.5*cm, 0, z_pos_lid), logicLid, "Phys_Lid", vacuum_stLV, false, 0, true);



  // ========================================================
  // EXTENSIÓN TORRE EJE Y+ (Horizontal)
  // ========================================================
  G4RotationMatrix* rotTowerY = new G4RotationMatrix();
  rotTowerY->rotateX(-90.*deg);

  G4double len_tower_y = half_h_shield; 
  G4double y_pos_tower = halfsize_box + len_tower_y; 


  // TUBO DE ACERO (EXTERIOR) ---
  /*
  G4Tubs* solid_Shield_Y = new G4Tubs("Shield_Y", RadioCyl_Tower - half_inch, RadioCyl_Tower, len_tower_y, 0, 2*pi);
  G4LogicalVolume* logicShieldY = new G4LogicalVolume(solid_Shield_Y, Steel, "Logic_Shield_Y");
  logicShieldY->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
  */
  new G4PVPlacement(rotTowerY, G4ThreeVector(0, y_pos_tower, 0), logicShieldSt, "Phys_Shield_St_Y", logicWorld, false, 0, true);


  // EL RELLENO DE VACÍO (INTERIOR) ---
  G4Tubs* solid_Vac_Y = new G4Tubs("Vac_Y", 0., RadioCyl_Tower - half_inch, len_tower_y, 0, 2*pi);
  G4LogicalVolume* logicVacY = new G4LogicalVolume(solid_Vac_Y, Vacuum, "Logic_Vac_Y");
  logicVacY->SetVisAttributes(G4VisAttributes(G4Colour(1., 1., 0., 0.3))); // Amarillo transparente
  new G4PVPlacement(rotTowerY, G4ThreeVector(0, y_pos_tower, 0), logicVacY, "Phys_Vac_Y", logicWorld, false, 0, true);


  // END CAP ---
  G4double y_pos_cap = y_pos_tower + len_tower_y + (half_inch/2);

  G4Tubs* solid_EndCapY = new G4Tubs("EndCapY", 0., RadioCyl_Tower, half_inch/2, 0., 2*pi);
  G4LogicalVolume* logicEndCapY = new G4LogicalVolume(solid_EndCapY, Steel, "Logic_EndCapY");
  logicEndCapY->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));

  new G4PVPlacement(rotTowerY, G4ThreeVector(0, y_pos_cap, 0), logicEndCapY, "Phys_EndCap_Y", logicWorld, false, 0, true);


  // ========================================================
  // STEEL FLAPS)
  // ========================================================

  G4double flap_thickness = half_inch / 2; // ~0.635 cm
  G4double flap_radius = halfsize_box;

  G4Tubs* solid_Flap = new G4Tubs("Flap", 0., flap_radius, flap_thickness, 0., 2*pi);
  G4LogicalVolume* logicFlap = new G4LogicalVolume(solid_Flap, Steel, "Logic_Flap");
  logicFlap->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));

  G4double dist_flap = halfsize_box + flap_thickness;

  
  // Flap +X
  new G4PVPlacement(rotY, G4ThreeVector(dist_flap, 0, 0), logicFlap, "Phys_Flap_PosX", logicWorld, false, 0, true);
  
  // Flap -X
  new G4PVPlacement(rotY, G4ThreeVector(-dist_flap, 0, 0), logicFlap, "Phys_Flap_NegX", logicWorld, false, 1, true);

  // Flap +Y
  // new G4PVPlacement(rotX, G4ThreeVector(0, dist_flap, 0), logicFlap, "Phys_Flap_PosY", logicWorld, false, 0, true);

  // Flap -Y
  new G4PVPlacement(rotX, G4ThreeVector(0, -dist_flap, 0), logicFlap, "Phys_Flap_NegY", logicWorld, false, 1, true);

  // Flap -Z
  new G4PVPlacement(0, G4ThreeVector(0, 0, -dist_flap), logicFlap, "Phys_Flap_NegZ", logicWorld, false, 0, true);

  // Flap +Z: there is not a flap here.

  // =======================================
  // END STEEL FLAPS
  // =======================================
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


