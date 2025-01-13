//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifdef G4ANALYSIS_USE
//#include "B02AnalysisManager.hh"
#endif

#include "B02SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "B02DetectorConstruction.hh"
#include "G4VUserTrackInformation.hh"
#include "G4VUserDetectorConstruction.hh"
#include "B02EventAction.hh"
#include "G4RunManager.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02SteppingAction::B02SteppingAction(B02EventAction* B02eventAction)
{
 fEventAction = B02eventAction;
}

void B02SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
 
       G4LogicalVolume *volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
      
      // Get the track associated with the step
       G4Track* track = aStep->GetTrack();
       auto name = track->GetParticleDefinition()->GetParticleName();
       auto trackID = track->GetTrackID(); 
       
       // get process, deposited energy and process name  
       const G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
       const G4VProcess* process = postStepPoint->GetProcessDefinedStep();
       auto edep = aStep->GetTotalEnergyDeposit();
       auto processName = process->GetProcessName();
       const B02DetectorConstruction *detectorConstruction = static_cast<const B02DetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
       
       G4LogicalVolume *fLogicalLAr = detectorConstruction->GetScoringVolume();
       // only give output in the active volume
       if(volume != fLogicalLAr)
        return;
    
    // step length associated with the primary muon 
    if (trackID == 1){
    G4double stepLength = aStep->GetStepLength();
    fEventAction->AddLength(stepLength);
    G4cout << "Track ID: " << trackID << " - Length traced on ActiveLAr: " << stepLength/cm << " cm" << G4endl; 
    G4cout << " particle type on active volume: " << name << G4endl; 
    G4cout << " Process name: " << processName << G4endl;
       //recolect the deposited energy of the muon that will decay 
         if (processName == "Decay"){
             fEventAction->AddEdep(edep);
             G4cout <<" Deposited Energy by Muon: " << edep/MeV << " MeV" << G4endl;
             }      
     }

}//steppingAction  
               
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

