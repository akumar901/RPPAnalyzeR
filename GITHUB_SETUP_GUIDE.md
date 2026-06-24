# How to put RPPAnalyzeR on GitHub — step by step

## Is this mine to publish? YES.

You commissioned this package for your data.
Your name (Amar Kumar) is listed as author in DESCRIPTION.
The MIT License lets you use, publish, and share it freely.
You can claim it fully as your own work.

---

## Part 1 — One-time computer setup (do this once ever)

### Install Git

**Windows:** Download from https://git-scm.com/download/win → install with defaults
**Mac:** Open Terminal → type `git --version` → if not installed, it will prompt you to install
**Linux:** `sudo apt install git`

### Configure Git with your name

Open Terminal (Mac/Linux) or Git Bash (Windows) and run:

```bash
git config --global user.name "Amar Kumar"
git config --global user.email "your.email@example.com"
```

### Create a GitHub account

Go to https://github.com → Sign Up → choose a username

---

## Part 2 — Create the GitHub repository

1. Log into GitHub
2. Click the **+** button (top right) → **New repository**
3. Fill in:
   - Repository name: `RPPAnalyzeR`
   - Description: `R pipeline for MD Anderson RPPA Core data analysis`
   - Set to **Public** (so others can install it with remotes::install_github)
   - Do NOT tick "Add a README" — we already have one
   - Do NOT tick "Add .gitignore" — we already have one
4. Click **Create repository**
5. GitHub will show you a page with setup instructions — keep it open

---

## Part 3 — Upload the package

Open Terminal (Mac/Linux) or Git Bash (Windows).

Navigate to the RPPAnalyzeR folder you downloaded:

```bash
cd /path/to/RPPAnalyzeR
```

(Replace `/path/to/` with wherever you unzipped the package, for example:
  Windows: `cd C:\Users\Amar\Downloads\RPPAnalyzeR`
  Mac:     `cd ~/Downloads/RPPAnalyzeR`)

Now run these commands one by one:

```bash
# 1. Initialise git in this folder
git init

# 2. Add all files
git add .

# 3. Make your first commit
git commit -m "Initial release: RPPAnalyzeR v0.1.0"

# 4. Name the branch 'main'
git branch -M main

# 5. Connect to GitHub (replace YOUR_USERNAME with your actual GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/RPPAnalyzeR.git

# 6. Push everything to GitHub
git push -u origin main
```

GitHub will ask for your username and password.
NOTE: GitHub no longer accepts your account password — use a Personal Access Token instead:
  - Go to GitHub → Settings → Developer Settings → Personal Access Tokens → Tokens (classic)
  - Click Generate new token → give it a name → tick "repo" → Generate token
  - Copy the token and use it as your password when Git asks

---

## Part 4 — Verify it worked

1. Go to `https://github.com/YOUR_USERNAME/RPPAnalyzeR`
2. You should see all your files and the README displayed automatically
3. The README.md will render beautifully with all the formatting

---

## Part 5 — Update your DESCRIPTION file

Open `DESCRIPTION` in any text editor and update:

```
Authors@R: 
    person("Amar", "Kumar", email = "your.real.email@example.com", role = c("aut", "cre"))
...
URL: https://github.com/YOUR_USERNAME/RPPAnalyzeR
BugReports: https://github.com/YOUR_USERNAME/RPPAnalyzeR/issues
```

Also update the badge URLs in README.md:
Replace `AmarKumar` with your actual GitHub username in the badge lines at the top.

Then push the update:

```bash
git add DESCRIPTION README.md
git commit -m "Update author email and GitHub URL"
git push
```

---

## Part 6 — Install your own package (test it works)

In R:

```r
remotes::install_github("YOUR_USERNAME/RPPAnalyzeR")
library(RPPAnalyzeR)

# Run your analysis
results <- run_pipeline(
  xlsx_path  = "01_Jonathan_Coloff__Vipin_Rawat.xlsx",
  output_dir = "my_rppa_output"
)
```

---

## Part 7 — Future updates (how to push changes)

Whenever you improve the package:

```bash
# See what changed
git status

# Add changed files
git add .

# Describe what you changed
git commit -m "Add pathway enrichment function"

# Push to GitHub
git push
```

---

## Useful GitHub features to enable

- **About section**: On your repo page, click the gear icon next to "About" → add a description and topics like `rppa`, `proteomics`, `r-package`, `bioinformatics`
- **Releases**: Go to Releases → Draft new release → tag it `v0.1.0` → this lets people install a specific version
- **Issues**: Leave issues enabled so collaborators can report bugs

---

## What to say in your publications

> "Data analysis was performed using the RPPAnalyzeR R package (Kumar A, 2025; https://github.com/YOUR_USERNAME/RPPAnalyzeR)."

And always include the MD Anderson RPPA Core grant acknowledgement:
> "The Functional Proteomics RPPA Core is supported by MD Anderson Cancer Center Support Grant # 5 P30 CA016672-40."
