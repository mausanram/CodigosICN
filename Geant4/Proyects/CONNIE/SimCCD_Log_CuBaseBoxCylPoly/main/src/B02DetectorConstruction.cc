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

  // Lead //
  auto Pb = nistManager->FindOrBuildMaterial("G4_Pb");

  // Polyethylene //
  auto Poly = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");

  // AlN //
  G4Element* elAl = nistManager->FindOrBuildElement("Al");
  G4Element* elN = nistManager->FindOrBuildElement("N");
  G4Material* AlN = new G4Material("AluminumNitride", 3.26 * g/cm3, 2);
  AlN->AddElement(elAl, 1);
  AlN->AddElement(elN, 1);

  // ====================== // end materials
    
  G4double units = 1 * cm; 
  // ==============================================
  // World
  // ===============================================
  fWorldLength = 200*units;
  G4double HalfWorldLength = 1.*fWorldLength;
  solidWorld= new G4Box("World",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  //logicWorld= new G4LogicalVolume( solidWorld, Vacuum, "World", 0, 0, 0);
  logicWorld= new G4LogicalVolume(solidWorld, Air, "World");
  physiWorld = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.), logicWorld, "World", 0, false, 0, true);   
  // =========================================== //
  
  // ============= Colors ================ //
  G4VisAttributes* attPb = new G4VisAttributes(G4Colour::Blue());
  // attSteel->SetForceAuxEdgeVisible(true);

  G4VisAttributes* attCu = new G4VisAttributes(G4Colour::Brown());
  // attCu->SetForceAuxEdgeVisible(true);
  // attCu->SetForceSolid(true);

  G4VisAttributes* attPoly = new G4VisAttributes(G4Colour::Cyan());
  // attAl->SetForceAuxEdgeVisible(true);
  // attAl->SetForceSolid(true);

  G4VisAttributes* attAlN = new G4VisAttributes(G4Colour::Gray());
  // attAlN->SetForceAuxEdgeVisible(true);
  // attAlN->SetForceSolid(true);

  G4VisAttributes* attSi = new G4VisAttributes(G4Colour::White());
  // attSi->SetForceAuxEdgeVisible(true);

  G4double fact_transparent = 0.1;
  G4VisAttributes* semiTransparentYellow = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, fact_transparent));

  fact_transparent = 0.7;
  G4VisAttributes* semiTransparentLighBlue = new G4VisAttributes(G4Colour(173., 216., 230., fact_transparent));

  G4VisAttributes* semiTransparentLighBrown = new G4VisAttributes(G4Colour(0.8, 0.6, 0.4, fact_transparent));
  semiTransparentLighBrown->SetForceAuxEdgeVisible(true);
  // semiTransparentLighBrown->SetForceSolid(true);
  // ====================================== //

  G4double half_inch = 2.54/2; // cm
  G4double tolerance = 0.0001 * units; 

  // ========= CU CYLINDER ======== //
  G4double cyl_high = 65/2 * units;
  G4double cyl_radius = 10. * units;

  G4double cyl_vac_high = cyl_high - half_inch * units;
  G4double cyl_vac_radius = cyl_radius - half_inch/2 * units;

  G4Tubs* solid_Cyl = new G4Tubs("solid_Cyl", 0., cyl_radius, cyl_high, 0., 2*pi);
  G4LogicalVolume* solid_CylCu = new G4LogicalVolume(solid_Cyl, Cu, "solid_CylCu");
  solid_CylCu->SetVisAttributes(semiTransparentLighBrown);

  G4Tubs* vac_Cyl = new G4Tubs("vac_Cyl", 0., cyl_vac_radius, cyl_vac_high, 0., 2*pi);
  G4LogicalVolume* vac_CylLog = new G4LogicalVolume(vac_Cyl, Vacuum, "vac_CylLog");
  vac_CylLog->SetVisAttributes(semiTransparentYellow);
  
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), solid_CylCu, "solid_CylCu", logicWorld, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(0, 0, cyl_high - cyl_vac_high), vac_CylLog, "vac_CylLog", solid_CylCu, false, 0, true);
  // ========= END CU CYLINDER ======== //

  // ==== CU FANGLE ==== //
  G4double flange_high = half_inch*1 * units;
  G4double flange_radius = cyl_radius + 48;

  G4double disk_high = flange_high;
  G4double disk_radius = flange_radius;

  G4double flange_vac_high = flange_high;
  G4double flange_vac_radius = cyl_vac_radius;

  G4Tubs* solid_flange = new G4Tubs("solid_flange", 0., flange_radius, flange_high, 0., 2*pi);
  G4LogicalVolume* logflangeCu = new G4LogicalVolume(solid_flange, Cu, "logflangeCu");
  logflangeCu->SetVisAttributes(semiTransparentLighBrown);

  G4Tubs* solid_disk = new G4Tubs("solid_disk", 0., disk_radius, disk_high, 0., 2*pi);
  G4LogicalVolume* logsolid_disk = new G4LogicalVolume(solid_disk, Cu, "logdisk");
  logsolid_disk->SetVisAttributes(semiTransparentLighBrown);

  G4Tubs* vac_flange= new G4Tubs("vac_Cyl", 0., flange_vac_radius, flange_vac_high, 0., 2*pi);
  G4LogicalVolume* logvac_flange= new G4LogicalVolume(vac_flange, Vacuum, "logvac_flange");
  logvac_flange->SetVisAttributes(semiTransparentYellow);

  new G4PVPlacement(0, G4ThreeVector(0, 0, cyl_high + flange_high), logflangeCu, "logflangeCu", logicWorld, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(0, 0, cyl_high + flange_high*3), logsolid_disk, "logsolid_disk", logicWorld, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logvac_flange, "logvac_flange", logflangeCu, false, 0, true);
  // ==== END CU FANGLE ==== //

  //==== LEAD CYLINDER ==== //
  G4double lead_cyl_high = 15/2 * units;
  G4double lead_cyl_radius = 9 * units;

  G4Tubs* solid_lead_cyl = new G4Tubs("solid_lead_cyl", 0., lead_cyl_radius, lead_cyl_high, 0., 2*pi);
  G4LogicalVolume* log_lead_cyl = new G4LogicalVolume(solid_lead_cyl, Pb, "log_lead_cyl");
  log_lead_cyl->SetVisAttributes(attPb);

  new G4PVPlacement(0, G4ThreeVector(0, 0, lead_cyl_high*5/2), log_lead_cyl, "log_lead_cyl", vac_CylLog, false, 0, true);


  // ========================================================
  // CU BOX
  // ========================================================
  G4double lateral_wall_x = 11.9/2 * units; // Units: cm
  G4double lateral_wall_y = 0.8/2 * units;
  G4double lateral_wall_z = 12.215/2 * units;

  G4double upper_wall_x = 11.9/2 * units;
  G4double upper_wall_y = 10.31/2 * units;
  G4double upper_wall_z = 0.9525/2 *units;

  G4double upper_base_x = 6.35/2 * units;
  G4double upper_base_y = 6.35/2 * units;
  G4double upper_base_z = 0.953/2 * units;

  G4Box* lateral_wall = new G4Box("lateral_wall", lateral_wall_x, lateral_wall_y, lateral_wall_z);
  G4LogicalVolume* lateral_wallCu = new G4LogicalVolume(lateral_wall, Cu, "lateral_wallCu");
  lateral_wallCu->SetVisAttributes(attCu);


  G4Box* upper_wall = new G4Box("upper_wall", upper_wall_x, upper_wall_y, upper_wall_z);
  G4LogicalVolume* upper_wallCu = new G4LogicalVolume(upper_wall, Cu, "upper_wallCu");
  upper_wallCu->SetVisAttributes(attCu);

  G4Box* upper_base = new G4Box("upper_base", upper_base_x, upper_base_y, upper_base_z);
  G4LogicalVolume* upper_baseCu = new G4LogicalVolume(upper_base, Cu, "upper_baseCu");
  upper_baseCu->SetVisAttributes(attCu);


  new G4PVPlacement(0, G4ThreeVector(0, -(upper_wall_y + lateral_wall_y + tolerance), 0), lateral_wallCu, "lateral_wallCu1", vac_CylLog, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(0, upper_wall_y+ lateral_wall_y+tolerance, 0), lateral_wallCu, "lateral_wallCu2", vac_CylLog, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(0, 0, lateral_wall_z - upper_wall_z + tolerance), upper_wallCu, "upper_wallCu", vac_CylLog, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(0, 0, lateral_wall_z + upper_base_z + tolerance), upper_baseCu, "upper_baseCu", vac_CylLog, false, 0, true);
  // ========================================================
  // END CU BOX
  // ========================================================

  // ========================================================
  // CU BASE  + CCD
  // ========================================================
  G4double z_offset = -(half_inch) * units; 

  /// Cu Base 
  G4double base_x = 11.1/2 * cm;
  G4double base_y = 10.2/2 * cm;
  G4double base_z = 0.2/2 * cm; // Grosor

  G4double cut_x = 7.722/2 * cm;
  G4double cut_y = 8.992/2 * cm;
  G4double cut_z = 0.089/2 * cm; // Profundidad corte

  // CCD
  G4double pixel_size = 0.0015 * cm; // 15 micras
  G4double ccd_x = 682 * pixel_size / 2;
  G4double ccd_y = 1022 * pixel_size / 2;
  G4double ccd_z = 0.0675 / 2 * cm;

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

  new G4PVPlacement(rotFlip, G4ThreeVector(0, 0, z_center_base), logicBaseCu, "Phys_BaseCu", vac_CylLog, false, 0, true);


  // Substrate (AlN)
  G4Box* solid_Substrate = new G4Box("Substrate", cut_x, cut_y, cut_z); // Mismo tamaño que el corte
  G4LogicalVolume* logicSubstrate = new G4LogicalVolume(solid_Substrate, AlN, "Logic_Substrate");
  logicSubstrate->SetVisAttributes(G4VisAttributes(G4Colour::Gray()));

  G4double z_substrate_pos = z_center_base - (base_z - cut_z); 
  new G4PVPlacement(rotFlip, G4ThreeVector(0, 0, z_substrate_pos), logicSubstrate, "Phys_Substrate", vac_CylLog, false, 0, true);


  G4Box* solid_CCD = new G4Box("CCD", ccd_x, ccd_y, ccd_z);
  G4LogicalVolume* logicCCD = new G4LogicalVolume(solid_CCD, Si, "Logic_CCD");
  logicCCD->SetVisAttributes(attSi); // Blanco brillante

  G4double z_ccd_pos = z_substrate_pos - cut_z - ccd_z;
  new G4PVPlacement(rotFlip, G4ThreeVector(0, 0, z_ccd_pos), logicCCD, "Phys_CCD", vac_CylLog, false, 0, true);


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

  // ===== POLYETHYLENE LAYER 1 ===== //
  G4double thickness_poly_layer = 30./2;
  G4int n_poly_walls = 10;
  G4double poly_thickness = (thickness_poly_layer/n_poly_walls) * units;

  G4double lateralwall_poly_x = poly_thickness;
  G4double lateralwall_poly_y = cyl_radius + (thickness_poly_layer*2)*units;
  G4double lateralwall_poly_z = 50/2 * units;

  // G4double lateralsmallwall_poly_x = lateralwall_poly_y - 60.1/2 * units;
  G4double lateralsmallwall_poly_x = cyl_radius;
  G4double lateralsmallwall_poly_y = poly_thickness;
  G4double lateralsmallwall_poly_z = lateralwall_poly_z;

  G4double lower_wall_poly_x = lateralwall_poly_y;
  G4double lower_wall_poly_y = lower_wall_poly_x;
  G4double lower_wall_poly_z = poly_thickness;

  G4Box* lateralwall_poly = new G4Box("lateralwall_poly", lateralwall_poly_x, lateralwall_poly_y, lateralwall_poly_z);
  G4LogicalVolume* logiclateralwall_poly = new G4LogicalVolume(lateralwall_poly, Poly, "logiclateralwall_poly");
  logiclateralwall_poly->SetVisAttributes(attPoly);

  G4Box* lateralsmallwall_poly = new G4Box("lateralsmallwall_poly", lateralsmallwall_poly_x, lateralsmallwall_poly_y, lateralsmallwall_poly_z);
  G4LogicalVolume* logiclateralsmallwall_poly = new G4LogicalVolume(lateralsmallwall_poly, Poly, "logiclateralsmallwall_poly");
  logiclateralsmallwall_poly->SetVisAttributes(attPoly);

  G4Box* lower_wall_poly = new G4Box("lower_wall_poly_x", lower_wall_poly_x, lower_wall_poly_y, lower_wall_poly_z);
  G4LogicalVolume* logiclower_wall_poly = new G4LogicalVolume(lower_wall_poly, Poly, "logiclower_wall_poly_x");
  logiclower_wall_poly->SetVisAttributes(attPoly);

  // Walls in X+ ---
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(cyl_radius + lateralwall_poly_x + 2*lateralwall_poly_x*i,0, -(cyl_high - lateralwall_poly_z)), logiclateralwall_poly, "logiclateralwall_poly", logicWorld, false, 0, true);  
  }
  // Walls in X- ---
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(-(cyl_radius + lateralwall_poly_x + 2*lateralwall_poly_x*i),0, -(cyl_high - lateralwall_poly_z)), logiclateralwall_poly, "logiclateralwall_poly", logicWorld, false, 0, true);  
  }

  // Walls in Y+ ---
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(0, cyl_radius + lateralsmallwall_poly_y + 2*lateralsmallwall_poly_y*i, -(cyl_high - lateralsmallwall_poly_z)), logiclateralsmallwall_poly, "logiclateralsmallwall_poly", logicWorld, false, 0, true);  
  }

  // Walls in Y- ---
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(0, -(cyl_radius + lateralsmallwall_poly_y + 2*lateralsmallwall_poly_y*i), -(cyl_high - lateralsmallwall_poly_z)), logiclateralsmallwall_poly, "logiclateralsmallwall_poly", logicWorld, false, 0, true);  
  }

  // Lower walls 
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(0,0, -(cyl_high + lower_wall_poly_z + 2*lower_wall_poly_z*i)), logiclower_wall_poly, "logiclower_wall_poly", logicWorld, false, 0, true);  
  }
  // ===== END POLYETHYLENE LAYER 1 ===== //
  
  // ===== LEAD LAYER ===== //
  G4int n_lead_brick = 3;
  G4double thickness_lead_layer = 15./2;
  G4double lead_thickness = (thickness_lead_layer/n_lead_brick) *units;

  G4double lateralwall_lead_x = lead_thickness;
  G4double lateralwall_lead_y = lateralwall_poly_y + thickness_poly_layer * units;
  G4double lateralwall_lead_z = cyl_high + thickness_poly_layer *units;

  G4double lateralsmallwall_lead_x = lateralwall_lead_y - thickness_poly_layer*units;
  G4double lateralsmallwall_lead_y = lead_thickness;
  G4double lateralsmallwall_lead_z = lateralwall_lead_z;

  G4double lowerwall_lead_x = lateralwall_lead_y;
  G4double lowerwall_lead_y = lateralwall_lead_y;
  G4double lowerwall_lead_z = lead_thickness;

  G4double upper_longwall_lead_x = lateralwall_lead_y - thickness_poly_layer*units;
  G4double upper_longwall_lead_y = 25.2/2*units;
  G4double upper_longwall_lead_z = lead_thickness;

  G4double upper_shortwall_lead_x = upper_longwall_lead_y;
  G4double upper_shortwall_lead_y = flange_radius;
  G4double upper_shortwall_lead_z = lead_thickness;

  G4Box* lateralwall_lead = new G4Box("lateralwall_lead", lateralwall_lead_x, lateralwall_lead_y, lateralwall_lead_z);
  G4LogicalVolume* logiclateralwall_lead = new G4LogicalVolume(lateralwall_lead, Pb, "logiclateralwall_lead");
  logiclateralwall_lead->SetVisAttributes(attPb);

  G4Box* lateralsmallwall_lead = new G4Box("lateralsmallwall_lead", lateralsmallwall_lead_x, lateralsmallwall_lead_y, lateralsmallwall_lead_z);
  G4LogicalVolume* logiclateralsmallwall_lead = new G4LogicalVolume(lateralsmallwall_lead, Pb, "logiclateralsmallwall_lead");
  logiclateralsmallwall_lead->SetVisAttributes(attPb);

  G4Box* lowerwall_lead = new G4Box("lowerwall_lead", lowerwall_lead_x, lowerwall_lead_y, lowerwall_lead_z);
  G4LogicalVolume* logiclowerwall_lead = new G4LogicalVolume(lowerwall_lead, Pb, "logiclowerwall_lead");
  logiclowerwall_lead->SetVisAttributes(attPb);

  G4Box* upper_longwall_lead = new G4Box("upper_longwall_lead", upper_longwall_lead_x, upper_longwall_lead_y, upper_longwall_lead_z);
  G4LogicalVolume* logicupper_longwall_lead = new G4LogicalVolume(upper_longwall_lead, Pb, "logicupper_longwall_lead");
  logicupper_longwall_lead->SetVisAttributes(attPb);

  G4Box* upper_shortwall_lead = new G4Box("upper_shortwall_lead", upper_shortwall_lead_x, upper_shortwall_lead_y, upper_shortwall_lead_z);
  G4LogicalVolume* logicupper_shortwall_lead = new G4LogicalVolume(upper_shortwall_lead, Pb, "logicupper_shortwall_lead");
  logicupper_shortwall_lead->SetVisAttributes(attPb);

  // // Walls in X+ ---
  for (int i = 0; i < n_lead_brick; i++){
    new G4PVPlacement(0, G4ThreeVector(cyl_radius + (thickness_poly_layer*2)*units + lateralwall_lead_x + 2*lateralwall_lead_x*i,0, -(thickness_poly_layer)*units), logiclateralwall_lead, "logiclateralwall_lead", logicWorld, false, 0, true);  
  }
  // Walls in X- ---
  for (int i = 0; i < n_lead_brick; i++){
    new G4PVPlacement(0, G4ThreeVector(-(cyl_radius + (thickness_poly_layer*2)*units + lateralwall_lead_x + 2*lateralwall_lead_x*i),0, -(thickness_poly_layer)*units), logiclateralwall_lead, "logiclateralwall_lead", logicWorld, false, 0, true);  
  }
  // Walls in Y+
  for (int i = 0; i < n_lead_brick; i++){
    new G4PVPlacement(0, G4ThreeVector(0, cyl_radius + (thickness_poly_layer*2)*units + lateralsmallwall_lead_y + 2*lateralsmallwall_lead_y*i, -(thickness_poly_layer)*units), logiclateralsmallwall_lead, "logiclateralsmallwall_lead", logicWorld, false, 0, true);  
  }
  // Walls in Y-
  for (int i = 0; i < n_lead_brick; i++){
    new G4PVPlacement(0, G4ThreeVector(0, -(cyl_radius + (thickness_poly_layer*2)*units + lateralsmallwall_lead_y + 2*lateralsmallwall_lead_y*i), -(thickness_poly_layer)*units), logiclateralsmallwall_lead, "logiclateralsmallwall_lead", logicWorld, false, 0, true);  
  }
  // Walls in Z-
  for (int i = 0; i < n_lead_brick; i++){
    new G4PVPlacement(0, G4ThreeVector(0, 0, -(cyl_high + (thickness_poly_layer*2)*units +  lowerwall_lead_z + 2*lowerwall_lead_z*i)), logiclowerwall_lead, "logiclowerwall_lead", logicWorld, false, 0, true);  
  }
  // Walls in Z+ long
  for (int i = 0; i < n_lead_brick; i++){
    new G4PVPlacement(0, G4ThreeVector(0, flange_radius + upper_longwall_lead_y, (lateralwall_poly_z - thickness_lead_layer*units + lowerwall_lead_z + 2*lowerwall_lead_z*i + 0.5*units)), logicupper_longwall_lead, "logicupper_longwall_lead", logicWorld, false, 0, true);  
  }
  // Walls in Z+ long
  for (int i = 0; i < n_lead_brick; i++){
    new G4PVPlacement(0, G4ThreeVector(0, -(flange_radius + upper_longwall_lead_y), (lateralwall_poly_z - thickness_lead_layer*units + lowerwall_lead_z + 2*lowerwall_lead_z*i + 0.5*units)), logicupper_longwall_lead, "logicupper_longwall_lead", logicWorld, false, 0, true);  
  }
  // Walls in Z+ short
  for (int i = 0; i < n_lead_brick; i++){
    new G4PVPlacement(0, G4ThreeVector((flange_radius + upper_shortwall_lead_x), 0, (lateralwall_poly_z - thickness_lead_layer*units + lowerwall_lead_z + 2*lowerwall_lead_z*i + 0.5*units)), logicupper_shortwall_lead, "logicupper_shortwall_lead", logicWorld, false, 0, true);  
  }
  // Walls in Z+ short
  for (int i = 0; i < n_lead_brick; i++){
    new G4PVPlacement(0, G4ThreeVector(-(flange_radius + upper_shortwall_lead_x - tolerance), 0, (lateralwall_poly_z - 15/2*units + lowerwall_lead_z + 2*lowerwall_lead_z*i)), logicupper_shortwall_lead, "logicupper_shortwall_lead", logicWorld, false, 0, true);  
  }
  // ===== END LEAD LAYER ===== //

  // ===== POLYETHYLENE LAYER 2 ===== //
  thickness_poly_layer = 30./2;
  n_poly_walls = 10;
  poly_thickness = (thickness_poly_layer/n_poly_walls) * units;

  lateralwall_poly_x = poly_thickness;
  lateralwall_poly_y = cyl_radius + (thickness_poly_layer*4 + thickness_lead_layer*2) * units;
  lateralwall_poly_z = cyl_high+ (thickness_poly_layer*2 + thickness_lead_layer)*units;

  lateralsmallwall_poly_x = lateralwall_poly_y - (thickness_poly_layer*2) * units;
  lateralsmallwall_poly_y = poly_thickness;
  lateralsmallwall_poly_z = lateralwall_poly_z;

  lower_wall_poly_x = lateralwall_poly_y;
  lower_wall_poly_y = lateralwall_poly_y;
  lower_wall_poly_z = poly_thickness;

  lateralwall_poly = new G4Box("lateralwall_lead", lateralwall_poly_x, lateralwall_poly_y, lateralwall_poly_z);
  logiclateralwall_poly = new G4LogicalVolume(lateralwall_poly, Poly, "logiclateralwall_poly");
  logiclateralwall_poly->SetVisAttributes(attPoly);

  lateralsmallwall_poly = new G4Box("lateralsmallwall_poly", lateralsmallwall_poly_x, lateralsmallwall_poly_y, lateralsmallwall_poly_z);
  logiclateralsmallwall_poly = new G4LogicalVolume(lateralsmallwall_poly, Poly, "logiclateralsmallwall_poly");
  logiclateralsmallwall_poly->SetVisAttributes(attPoly);

  lower_wall_poly = new G4Box("lower_wall_poly_x", lower_wall_poly_x, lower_wall_poly_y, lower_wall_poly_z);
  logiclower_wall_poly = new G4LogicalVolume(lower_wall_poly, Poly, "logiclower_wall_poly_x");
  logiclower_wall_poly->SetVisAttributes(attPoly);

  // Walls in X+ ---
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(cyl_radius + (thickness_poly_layer*2 + thickness_lead_layer*2 +0.1)*units +lateralwall_poly_x + 2*lateralwall_poly_x*i,0, -(cyl_high- (thickness_poly_layer + thickness_lead_layer + half_inch + half_inch/2)*units)), logiclateralwall_poly, "logiclateralwall_poly", logicWorld, false, 0, true);  
  }
  // Walls in X- ---
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(-(cyl_radius + (thickness_poly_layer*2 + thickness_lead_layer*2 +0.1)*units +lateralwall_poly_x + 2*lateralwall_poly_x*i),0, -(cyl_high- (thickness_poly_layer + thickness_lead_layer + half_inch + half_inch/2)*units)), logiclateralwall_poly, "logiclateralwall_poly", logicWorld, false, 0, true);  
  }
  // Walls in Y+ ---
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(0, cyl_radius + (thickness_poly_layer*2 + thickness_lead_layer*2 +0.1)*units +lateralsmallwall_poly_y + 2*lateralsmallwall_poly_y*i, -(cyl_high- (thickness_poly_layer + thickness_lead_layer + half_inch + half_inch/2)*units)), logiclateralsmallwall_poly, "logiclateralsmallwall_poly", logicWorld, false, 0, true);  
  }
  // Walls in Y- ---
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(0, -(cyl_radius + (thickness_poly_layer*2 + thickness_lead_layer*2 +0.1)*units +lateralsmallwall_poly_y + 2*lateralsmallwall_poly_y*i), -(cyl_high- (thickness_poly_layer + thickness_lead_layer + half_inch + half_inch/2)*units)), logiclateralsmallwall_poly, "logiclateralsmallwall_poly", logicWorld, false, 0, true);  
  }
  // Lower walls 
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(0,0, -(cyl_high + (thickness_poly_layer*2 + thickness_lead_layer*2 +0.1)*units + lower_wall_poly_z + 2*lower_wall_poly_z*i)), logiclower_wall_poly, "logiclower_wall_poly", logicWorld, false, 0, true);  
  }


  // Upper walls 
  lateralwall_poly_x = poly_thickness;
  lateralwall_poly_y = flange_radius + (thickness_poly_layer*2 + 20.5/2) * units;
  lateralwall_poly_z = cyl_high+ (0)*units;

  lateralsmallwall_poly_x = lateralwall_poly_y - (thickness_poly_layer*2) * units;
  lateralsmallwall_poly_y = poly_thickness;
  lateralsmallwall_poly_z = lateralwall_poly_z;

  G4double upper_wall_poly_x = lateralwall_poly_y;
  G4double upper_wall_poly_y = lateralwall_poly_y;
  G4double upper_wall_poly_z = poly_thickness;

  lateralwall_poly = new G4Box("lateralwall_lead", lateralwall_poly_x, lateralwall_poly_y, lateralwall_poly_z);
  logiclateralwall_poly = new G4LogicalVolume(lateralwall_poly, Poly, "logiclateralwall_poly");
  logiclateralwall_poly->SetVisAttributes(attPoly);

  lateralsmallwall_poly = new G4Box("lateralsmallwall_poly", lateralsmallwall_poly_x, lateralsmallwall_poly_y, lateralsmallwall_poly_z);
  logiclateralsmallwall_poly = new G4LogicalVolume(lateralsmallwall_poly, Poly, "logiclateralsmallwall_poly");
  logiclateralsmallwall_poly->SetVisAttributes(attPoly);

  G4Box* upperwall_poly = new G4Box("upperwall_poly", upper_wall_poly_x, upper_wall_poly_y, upper_wall_poly_z);
  G4LogicalVolume*  logicupperwall_poly = new G4LogicalVolume(upperwall_poly, Poly, "logicupperwall_poly");
  logicupperwall_poly->SetVisAttributes(attPoly);

  // Walls in X+
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector( flange_radius + (20.5/2) * units+lateralwall_poly_x + 2*lateralwall_poly_x*i,0, (cyl_high*2+flange_high*4 + (1.)*units)), logiclateralwall_poly, "logiclateralwall_poly", logicWorld, false, 0, true);  
  }
  // Walls in X-
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector( -(flange_radius + (20.5/2) * units+lateralwall_poly_x + 2*lateralwall_poly_x*i),0, (cyl_high*2+flange_high*4 + (1.)*units)), logiclateralwall_poly, "logiclateralwall_poly", logicWorld, false, 0, true);  
  }
  // Walls in Y+
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(0, flange_radius + (20.5/2) * units +lateralsmallwall_poly_y + 2*lateralsmallwall_poly_y*i, (cyl_high*2+flange_high*4 + (1.)*units)), logiclateralsmallwall_poly, "logiclateralsmallwall_poly", logicWorld, false, 0, true);  
  }
  // Walls in Y- 
  for (int i = 0; i < n_poly_walls; i++){
    new G4PVPlacement(0, G4ThreeVector(0, -(flange_radius + (20.5/2) * units +lateralsmallwall_poly_y + 2*lateralsmallwall_poly_y*i), (cyl_high*2+flange_high*4 + (1.)*units)), logiclateralsmallwall_poly, "logiclateralsmallwall_poly", logicWorld, false, 0, true);  
  }

  new G4PVPlacement(0, G4ThreeVector(0, 0, (cyl_high*2 + flange_high*4 + lateralwall_poly_z*1 + upper_wall_poly_z*2)), logicupperwall_poly, "logicupperwall_poly", logicWorld, false, 0, true);  

  // ===== END POLYETHYLENE LAYER 2 ===== //

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
  // new G4PVPlacement(rotFlip, G4ThreeVector(x_pos_conn1, 0, z_pos_conn1), logicConn1, "Phys_Conn1", logicWorld, false, 0, true);


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

  // new G4PVPlacement(rotFlip, G4ThreeVector(x_pos_conn2, 0, z_pos_conn2), logicConn2, "Phys_Conn2", logicWorld, false, 0, true);


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
  // new G4PVPlacement(rotFlip, G4ThreeVector(-0.5*cm, 0, z_pos_lid), logicLid, "Phys_Lid", logicWorld, false, 0, true);
  
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


