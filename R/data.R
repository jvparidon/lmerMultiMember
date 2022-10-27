#' National Football League scores from the 2021 season
#'
#' A subset of NFL scores from the FiveThirtyEight's football model dataset,
#' released under CC-BY license and reproduced here with attribution.
#'
#' @format ## `nfl_scores_2021`
#' A data frame with 285 rows and 7 columns:
#' \describe{
#'   \item{date}{Date on which a game was played}
#'   \item{season}{NFL season in which a game was played}
#'   \item{home_team, visiting_team}{3-letter codes for the playing teams}
#'   \item{home_score, visiting_score}{Scores for the playing teams}
#'   \item{winner}{Whether the home team or the visiting team won the game}
#' }
#' @source <https://github.com/fivethirtyeight/data/tree/master/nfl-elo>
"nfl_scores_2021"
