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

  // AlN //
  G4Element* elAl = nistManager->FindOrBuildElement("Al");
  G4Element* elN = nistManager->FindOrBuildElement("N");
  G4Material* AlN = new G4Material("AluminumNitride", 3.26 * g/cm3, 2);
  AlN->AddElement(elAl, 1);
  AlN->AddElement(elN, 1);
  // ====================== // end materials
    
  // ============== WORLD ============== //
  fWorldLength = 50*cm;
  G4double HalfWorldLength = 1.*fWorldLength;
  solidWorld= new G4Box("World",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  //logicWorld= new G4LogicalVolume( solidWorld, Vacuum, "World", 0, 0, 0);
  logicWorld= new G4LogicalVolume(solidWorld, Air, "World");
  physiWorld = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.), logicWorld, "World", 0, false, 0, true);   
  // =========================================== //
  
  // ============= Colors ================ //
  G4VisAttributes cgray(G4Colour::Gray());
  G4VisAttributes brown(G4Colour::Brown());
  G4VisAttributes white(G4Colour::White());

  G4VisAttributes* attSteel = new G4VisAttributes(G4Colour::Blue());
  attSteel->SetForceAuxEdgeVisible(true);

  G4double fact_transparent = 0.9;
  G4VisAttributes* semiTransparentYellow = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, fact_transparent));
  semiTransparentYellow->SetForceAuxEdgeVisible(true);
  semiTransparentYellow->SetForceSolid(true);

  G4VisAttributes* semiTransparentLighBlue = new G4VisAttributes(G4Colour(173., 216., 230., fact_transparent));
  // ================================= //


  // ============= STEEL BOX ============== //
  G4double half_inch = 2.54/2; // cm
  G4double halfsize_box = 12.7 * cm;      
  G4double thick_wall = half_inch;
  G4double inner_halfsize = halfsize_box - thick_wall; 
  G4double RadioCyl = halfsize_box - half_inch; 

  G4Box* solid_Box_Outer = new G4Box("Box_Outer", halfsize_box, halfsize_box, halfsize_box);
  G4Box* solid_Box_Inner = new G4Box("Box_Inner", inner_halfsize, inner_halfsize, inner_halfsize);

  G4SubtractionSolid* solid_FinalSteel = new G4SubtractionSolid("Steel_Final", solid_Box_Outer, solid_Box_Inner);

  G4LogicalVolume* logicSteelBox = new G4LogicalVolume(solid_FinalSteel, Steel, "Logic_SteelBox");
  new G4PVPlacement(0, G4ThreeVector(0,0,0), logicSteelBox, "Phys_SteelBox", logicWorld, false, 0, true);
  logicSteelBox->SetVisAttributes(attSteel);

  G4LogicalVolume* logicVacuumInterior = new G4LogicalVolume(solid_Box_Inner, Vacuum, "Logic_VacuumInterior");
  new G4PVPlacement(0, G4ThreeVector(0,0,0), logicVacuumInterior, "Phys_VacuumInterior", logicWorld, false, 0, true);
  logicVacuumInterior->SetVisAttributes(semiTransparentYellow);

  G4LogicalVolume* vacuum_stLV = logicVacuumInterior;
  // ============== END STEEL BOX ========== //

  // =================== CU BASE  + CCD =====================
  G4double z_offset = 0. * cm; 

  /// Cu Base 
  G4double base_x = 11.1/2 * cm;
  G4double base_y = 10.2/2 * cm;
  G4double base_z = 0.2/2 * cm; // Grosor

  G4double cut_x = 7.722/2 * cm;
  G4double cut_y = 8.992/2 * cm;
  G4double cut_z = 0.089/2 * cm; // Profundidad corte

  // CCD
  G4double pixel_size = 0.0015 * cm; // 15 micras
  G4double ccd_x = 250 * pixel_size / 2;
  G4double ccd_y = 529 * pixel_size / 2;
  G4double ccd_z = 0.0725 / 2 * cm;

  G4double z_start_base = (base_z - cut_z) + cut_z + ccd_z + base_z + z_offset;

  G4RotationMatrix* rotFlip = new G4RotationMatrix();
  rotFlip->rotateX(180.*deg);

  // Base 1
  G4Box* solid_BaseMain = new G4Box("BaseMain", base_x, base_y, base_z);

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

  // ----------- DETALLES INTERNOS (CONECTORES Y TAPA)

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
  // =================== END CU BASE  + CCD =====================
  
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


