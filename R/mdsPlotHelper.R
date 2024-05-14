calcSquareControlPoints <- function (x1, y1, x2, y2, curvature, angle, ncp) {
  dx <- x2 - x1
  dy <- y2 - y1
  slope <- dy/dx
  end <- (slope > 1 | (slope < 0 & slope > -1))
  if (curvature < 0) 
    end <- !end
  startx <- ifelse(end, x1, ifelse(abs(slope) > 1, x2 - dx, 
                                   x2 - sign(slope) * dy))
  starty <- ifelse(end, y1, ifelse(abs(slope) > 1, y2 - sign(slope) * 
                                     dx, y2 - dy))
  endx <- ifelse(end, ifelse(abs(slope) > 1, x1 + dx, x1 + 
                               sign(slope) * dy), x2)
  endy <- ifelse(end, ifelse(abs(slope) > 1, y1 + sign(slope) * 
                               dx, y1 + dy), y2)
  cps <- calcControlPoints(startx, starty, endx, endy, curvature, 
                           angle, ncp)
  ncurve <- length(x1)
  cps$x <- interleave(ncp, ncurve, cps$x, startx, endx, end)
  cps$y <- interleave(ncp, ncurve, cps$y, starty, endy, end)
  list(x = cps$x, y = cps$y, end = end)
}


calcControlPoints <- function (x1, y1, x2, y2, curvature, angle, ncp) {
  xm <- (x1 + x2)/2
  ym <- (y1 + y2)/2
  dx <- x2 - x1
  dy <- y2 - y1
  slope <- dy/dx
  if (is.null(angle)) {
    angle <- ifelse(slope < 0, 2 * atan(abs(slope)), 2 * 
                      atan(1/slope))
  }
  else {
    angle <- angle/180 * pi
  }
  sina <- sin(angle)
  cosa <- cos(angle)
  cornerx <- xm + (x1 - xm) * cosa - (y1 - ym) * sina
  cornery <- ym + (y1 - ym) * cosa + (x1 - xm) * sina
  beta <- -atan((cornery - y1)/(cornerx - x1))
  sinb <- sin(beta)
  cosb <- cos(beta)
  newx2 <- x1 + dx * cosb - dy * sinb
  newy2 <- y1 + dy * cosb + dx * sinb
  scalex <- (newy2 - y1)/(newx2 - x1)
  newx1 <- x1 * scalex
  newx2 <- newx2 * scalex
  ratio <- 2 * (sin(atan(curvature))^2)
  origin <- curvature - curvature/ratio
  if (curvature > 0) 
    hand <- "right"
  else hand <- "left"
  oxy <- calcOrigin(newx1, y1, newx2, newy2, origin, hand)
  ox <- oxy$x
  oy <- oxy$y
  dir <- switch(hand, left = -1, right = 1)
  maxtheta <- pi + sign(origin * dir) * 2 * atan(abs(origin))
  theta <- seq(0, dir * maxtheta, dir * maxtheta/(ncp + 1))[c(-1, 
                                                              -(ncp + 2))]
  costheta <- cos(theta)
  sintheta <- sin(theta)
  cpx <- ox + ((newx1 - ox) %*% t(costheta)) - ((y1 - oy) %*% 
                                                  t(sintheta))
  cpy <- oy + ((y1 - oy) %*% t(costheta)) + ((newx1 - ox) %*% 
                                               t(sintheta))
  cpx <- cpx/scalex
  sinnb <- sin(-beta)
  cosnb <- cos(-beta)
  finalcpx <- x1 + (cpx - x1) * cosnb - (cpy - y1) * sinnb
  finalcpy <- y1 + (cpy - y1) * cosnb + (cpx - x1) * sinnb
  list(x = as.numeric(t(finalcpx)), y = as.numeric(t(finalcpy)))
}

calcOrigin <- function (x1, y1, x2, y2, origin, hand) {
  xm <- (x1 + x2)/2
  ym <- (y1 + y2)/2
  dx <- x2 - x1
  dy <- y2 - y1
  slope <- dy/dx
  oslope <- -1/slope
  tmpox <- ifelse(!is.finite(slope), xm, ifelse(!is.finite(oslope), 
                                                xm + origin * (x2 - x1)/2, xm + origin * (x2 - x1)/2))
  tmpoy <- ifelse(!is.finite(slope), ym + origin * (y2 - y1)/2, 
                  ifelse(!is.finite(oslope), ym, ym + origin * (y2 - y1)/2))
  sintheta <- -1
  ox <- xm - (tmpoy - ym) * sintheta
  oy <- ym + (tmpox - xm) * sintheta
  list(x = ox, y = oy)
}

interleave <- function (ncp, ncurve, val, sval, eval, e) {
  sval <- rep(sval, length.out = ncurve)
  eval <- rep(eval, length.out = ncurve)
  result <- matrix(NA, ncol = ncurve, nrow = ncp + 1)
  m <- matrix(val, ncol = ncurve)
  for (i in 1L:ncurve) {
    if (e[i]) 
      result[, i] <- c(m[, i], eval[i])
    else result[, i] <- c(sval[i], m[, i])
  }
  as.numeric(result)
}
