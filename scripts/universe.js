// const baseURL = "http://127.0.0.1:3000"; // local
const baseURL = "https://qte7wuo072.execute-api.eu-west-2.amazonaws.com/Prod"; // cloud

const image = () => {
  // Definitions
  const url = baseURL + "/nbody";

  // Cosmological parameters
  const Omega_m = document.getElementById("Omega_m").value;
  const Omega_b = document.getElementById("Omega_b").value;
  const H_0 = document.getElementById("H_0").value;
  const sigma_8 = 0.8;
  // const sigma_8 = document.getElementById("sigma_8").value;
  const A_s = 2e-9;
  const n_s = document.getElementById("n_s").value;
  // const w_0 = -1;
  const w_0 = document.getElementById("w_0").value;
  const w_a = 0;
  const m_nu = 0;

  // Color scheme
  const color = "digilab";

  // Constants
  const kmin = 1e-3;
  const kmax = 1e1;
  const nk = 100;
  const z = 0;
  const npix = 1000;
  const Lbox = 500;
  const Tbox = 1;

  // Construct json
  const params = {
    method: "POST", // Unless this is present it will default to "GET"
    headers: {
      "Content-Type": "application/json",
      Accept: "application/json",
    },
    body: JSON.stringify({
      kmin: kmin,
      kmax: kmax,
      nk: nk,
      z: z,
      color: color,
      npix: npix,
      Lbox: Lbox,
      Tbox: Tbox,
      Omega_m: Omega_m,
      Omega_b: Omega_b,
      H_0: H_0,
      sigma_8: sigma_8,
      A_s: A_s,
      n_s: n_s,
      w_0: w_0,
      w_a: w_a,
      m_nu: m_nu,
    }),
  };

  // Show spinner and grey overlay
  const spinner = document.getElementById("spinner");
  const overlay = document.getElementById("overlay");
  const button = document.getElementById("clickButton");
  spinner.style.display = "block";
  overlay.style.display = "block";
  button.disabled = true;

  // Fetch
  console.log("Request sent");
  fetch(url, params)
    .then((response) => response.json()) // TODO: Not blob; maybe use response.blob()?
    .then((blob) => {
      console.log("Response blob received");
      console.log("Blob: " + blob);
      const data = blob.data;
      const image = "data:image/png;base64," + data;
      document.getElementById("image").src = image; // To set image within html
      console.log("Image displayed");
      spinner.style.display = "none";
      overlay.style.display = "none";
      button.disabled = false;
    })
    .catch((error) => {
      console.log("Error:", error);
      console.log("Failed to sample image");
      spinner.style.display = "none";
      overlay.style.display = "none";
      button.disabled = false;
    });
};
