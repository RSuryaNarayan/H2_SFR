#include "prob.H"

std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"') {
      count++;
    }
  }
  return (count % 2) == 0;
}

void
read_pmf(const std::string& myfile)
{
  std::string firstline;
  std::string secondline;
  std::string remaininglines;
  unsigned int pos1;
  unsigned int pos2;
  int variable_count;
  int line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  pos1 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    pos1 = pos2 + 1;
  }

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  // for (int i = 0; i < variable_count; i++)
  //  amrex::Print() << "Variable found: " << pmf_names[i] <<
  //  std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  PeleC::h_prob_parm_device->pmf_N = line_count;
  PeleC::h_prob_parm_device->pmf_M = variable_count - 1;
  PeleC::prob_parm_host->h_pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->h_pmf_Y.resize(
    static_cast<long>(PeleC::h_prob_parm_device->pmf_N) *
    PeleC::h_prob_parm_device->pmf_M);
  PeleC::prob_parm_host->pmf_Y.resize(
    static_cast<long>(PeleC::h_prob_parm_device->pmf_N) *
    PeleC::h_prob_parm_device->pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (int i = 0; i < PeleC::h_prob_parm_device->pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PeleC::prob_parm_host->h_pmf_X[i];
    for (int j = 0; j < PeleC::h_prob_parm_device->pmf_M; j++) {
      sinput >> PeleC::prob_parm_host
                  ->h_pmf_Y[j * PeleC::h_prob_parm_device->pmf_N + i];
    }
  }

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_X.begin(),
    PeleC::prob_parm_host->h_pmf_X.end(), PeleC::prob_parm_host->pmf_X.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_Y.begin(),
    PeleC::prob_parm_host->h_pmf_Y.end(), PeleC::prob_parm_host->pmf_Y.begin());
  PeleC::h_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::h_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
}


void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
 // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("do_perturb", PeleC::h_prob_parm_device->do_perturb);
    pp.query("P_amb", PeleC::h_prob_parm_device->P_amb);
    pp.query("T_amb", PeleC::h_prob_parm_device->T_amb);
    pp.query("Y_H2", PeleC::h_prob_parm_device->Y_H2);
    pp.query("Y_O2", PeleC::h_prob_parm_device->Y_O2);
    pp.query("Y_N2", PeleC::h_prob_parm_device->Y_N2);
    pp.query("D_CJ", PeleC::h_prob_parm_device->D_CJ);
    pp.query("frac_rho_fluc", PeleC::h_prob_parm_device->frac_rho_fluc);
    pp.query("delta_half", PeleC::h_prob_parm_device->delta_half);
    pp.query("PostStep", PeleC::h_prob_parm_device->PostStep);
    pp.query("x_ZND", PeleC::h_prob_parm_device->x_ZND);
    pp.query("A", PeleC::h_prob_parm_device->A);
    pp.query("n", PeleC::h_prob_parm_device->n);
    pp.query("t_end_turb", PeleC::h_prob_parm_device->t_end_turb);
    
    std::string znd_datafile;
    pp.query("znd_datafile", znd_datafile);
    read_pmf(znd_datafile);
  }
  PeleC::h_prob_parm_device->massfrac[H2_ID] = PeleC::h_prob_parm_device->Y_H2;
  PeleC::h_prob_parm_device->massfrac[O2_ID] = PeleC::h_prob_parm_device->Y_O2;
  PeleC::h_prob_parm_device->massfrac[N2_ID] = PeleC::h_prob_parm_device->Y_N2;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(PeleC::h_prob_parm_device->P_amb,
    	       PeleC::h_prob_parm_device->massfrac.begin(),
               PeleC::h_prob_parm_device->T_amb, 
               PeleC::h_prob_parm_device->rho_amb,
               PeleC::h_prob_parm_device->eint_amb);
  amrex::Print()<<"\nRHO QUIESCENT= "<<PeleC::h_prob_parm_device->rho_amb;
  amrex::Print()<<"\nINTERNAL ENERGY QUIESCENT= "<<PeleC::h_prob_parm_device->eint_amb<<"\n";
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
