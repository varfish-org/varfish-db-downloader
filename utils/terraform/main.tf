# Mangement of the GitHub project.

resource "github_repositoryvarfish-db-downloader" {
  name        = "varfish-db-downloader"
  description = "Download public databases for VarFish"

  has_issues = true
  visibility = "public"

  allow_rebase_merge     = false
  allow_merge_commit     = false
  delete_branch_on_merge = true

  vulnerability_alerts = true

  squash_merge_commit_message = "BLANK"
  squash_merge_commit_title   = "PR_TITLE"
}
